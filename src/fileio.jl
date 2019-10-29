# import Base: parse

include("qe/fileio.jl")
include("abinit/fileio.jl")
include("wannier90/fileio.jl")
include("elk/fileio.jl")

#--------------------Used by other file processing------------------#
function parse_k_line(line, T)
    splt = split(line)
    k1   = parse(T, splt[5])
    k2   = parse(T, splt[6])
    k3   = parse(T, splt[7][1:1:end-2])
    return Vec3([k1, k2, k3])
end

function write_flag_line(f, flag, data, seperator="=", i="")
    flagstr = string(flag)
    if flagstr[end-1] == '_' && tryparse(Int, string(flagstr[end])) != nothing
        flagstr = flagstr[1:end-2] * "($(flagstr[end]))"
    end
    write(f,"  $flagstr$i $seperator ")

    if typeof(data) <: Array

        if length(data) < 20
            write(f, "  $(data[1])")
            for x in data[2:end]
                write(f, " $x")
            end
            write(f, "\n")
        else
            write(f, "\n")
            for i = 1:3:length(data)
                write(f, "  $(data[i]) $(data[i + 1]) $(data[i + 2])\n")
            end
        end

    else #this should work for anything singular valued data such as bools, ''s and other types
        write(f, "$data\n")
    end

end

function parse_flag_val(val, T=Float64)
    if T == String
        return val
    end
    if occursin("d", val) && T <: Number
        val = replace(val, "d" => "e")
    end

    val = strip(val, '.')
    t = parse.(eltype(T), split(lowercase(val)))
    #deal with abinit constants -> all flags that are read which are not part of the abi[:structure] get cast into the correct atomic units!
    if length(t) > 1 && typeof(t[end]) == Symbol
        t = t[1:end-1] .* abi_conversions[t[end]]
    end

    length(t) == 1 ? t[1] : typeof(t) <: Vector{Real} ? convert.(T,t) : t
end

function write_data(f, data)
    if typeof(data) <: Vector{Vector{Float64}} || typeof(data) <: Vector{NTuple{4, Float64}} #k_points
        for x in data
            for y in x
                write(f, " $y")
            end
            write(f, "\n")
        end
    elseif typeof(data) <: Vector{Int} || typeof(data) <: NTuple{6, Int}
        for x in data
            write(f, " $x")
        end
        write(f, "\n")
    elseif typeof(data) <: Matrix
        im, jm = size(data)
        for i = 1:im
            for j = 1:jm
                write(f, " $(data[i, j])")
            end
            write(f, "\n")
        end
    end
end

Base.length(::Type{<:Number}) = 1

function Base.parse(::Type{NamedTuple{names, types}}, spl::Vector{<:AbstractString}) where {names, types}
    @assert sum(length.(types.parameters)) == length(spl)
    tbuf = Vector{Union{types.parameters...}}(undef, length(types.parameters))
    counter = 1
    for (i, typ) in enumerate(types.parameters)
        if typ <: AbstractVector
            l = length(typ)
            tbuf[i] = parse(typ, spl[counter:counter+l-1])
            counter += l
        else
            tbuf[i] = parse(typ, spl[counter])
            counter += 1
        end
    end
    pstring = (tbuf...,)
    return NamedTuple{names}(pstring)
end
Base.parse(::Type{T}, s::AbstractString) where {T <: NamedTuple} = parse(T, split(s))

function Base.parse(::Type{Point{N, T}}, spl::Vector{<:AbstractString}) where {N, T}
    @assert N == length(spl)
    return Point{N, T}(parse.(T, spl))
end
Base.parse(::Type{T}, s::AbstractString) where {T <: Point} = parse(T, split(s))

#---------------------------BEGINNING GENERAL SECTION-------------------#
#Incomplete: only works with SBATCH right now
function write_job_name(f, job::DFJob)
    write(f, "#SBATCH -J $(job.name) \n")
end

function write_job_header(f, job::DFJob)
    job_header = job.header == "" ? getdefault_jobheader() : job.header
    for line in job_header
        if occursin("\n", line)
            write(f, line)
        else
            write(f, line * "\n")
        end
    end
end

function writetojob(f, job, inputs::Vector{DFInput{Abinit}}; kwargs...)
    abinit_jobfiles   = write_abi_datasets(inputs, job.local_dir; kwargs...)
    abifiles = String[]
    num_abi = 0
    for (filename, pseudos, runcommand) in abinit_jobfiles
        push!(abifiles, filename)
        file, ext = splitext(filename)
        write(f, "$runcommand << !EOF\n$filename\n$(file * ".out")\n$(job.name * "_Xi$num_abi")\n$(job.name * "_Xo$num_abi")\n$(job.name * "_Xx$num_abi")\n")
        for pp in pseudos
            write(f, "$pp\n")
        end
        write(f, "!EOF\n")
        num_abi += 1
    end
    return abifiles
end

function writetojob(f, job, inputs::Vector{DFInput{Elk}}; kwargs...)
    save(inputs, job.structure; kwargs...)
    should_run = any(map(x->x.run, inputs)) 
    if !should_run
        write(f, "#")
    end
    writeexec.((f,), execs(inputs[1]))
    write(f, "< $(infilename(inputs[1])) > $(outfilename(inputs[1]))\n")
    return inputs
end

function writeexec(f, exec::Exec)
    direxec = joinpath(exec.dir, exec.exec)
    write(f, "$direxec")
    for flag in exec.flags
        write(f, " -$(flag.symbol)")
        if !isa(flag.value, AbstractString)
            for v in flag.value
                write(f," $v")
            end
        else
            write(f, " $(flag.value)")
        end
    end
    write(f, " ")
end

function writetojob(f, job, input::DFInput; kwargs...)
    filename    = infilename(input)
    should_run  = input.run
    save(input, job.structure; kwargs...)
    if !should_run
        write(f, "#")
    end
    writeexec.((f, ), execs(input))
    write(f, "< $filename > $(outfilename(input))\n")
    return (input,)
end

function writetojob(f, job, _input::DFInput{Wannier90}; kwargs...)
    filename    = infilename(_input)
    should_run  = _input.run
    id = findfirst(isequal(_input), job.inputs)
    seedname = name(_input)

	nscf_calc = getfirst(x -> isnscfcalc(x), job.inputs)
    runexec   = nscf_calc.execs
    # For elk the setup necessary for the wan_calc needs to be done before writing the wan input
    # because it's inside elk.in
	if package(nscf_calc) == QE
	    pw2waninput = qe_generate_pw2waninput(_input, "$(job.name)", runexec)
	    preprocess  = pop!(flags(_input), :preprocess)

	    if !preprocess || !should_run
	        write(f, "#")
	    end
	    writeexec.((f,), execs(_input))
	    write(f, "-pp $filename > $(outfilename(_input))\n")

	    save(_input, job.structure; kwargs...)
	    writetojob(f, job, pw2waninput; kwargs...)
	    flags(_input)[:preprocess] = preprocess
    elseif package(nscf_calc) == Elk
	    pw2waninput = job["elk2wannier"]
    end

    if !should_run
        write(f, "#")
    end
    writeexec.((f, ), execs(_input))
    write(f, "$filename > $(outfilename(_input))\n")
    return _input, pw2waninput
end

"""
    writejobfiles(job::DFJob; kwargs...)

Writes all the input files and job file that are linked to a DFJob.
Kwargs will be passed down to various writetojob functions.
"""
function writejobfiles(job::DFJob; kwargs...)
    rm.(joinpath.(Ref(job.local_dir), searchdir(job.local_dir, ".in")))
    open(joinpath(job.local_dir, "job.tt"), "w") do f
        write(f, "#!/bin/bash\n")
        write_job_name(f, job)
        write_job_header(f, job)
        written_inputs = DFInput[]
        abiinputs = Vector{DFInput{Abinit}}(filter(x -> package(x) == Abinit, inputs(job)))
        !isempty(abiinputs) && writetojob(f, job, abiinputs; kwargs...)
        elkinputs = Vector{DFInput{Elk}}(filter(x -> package(x) == Elk, inputs(job)))
        !isempty(elkinputs) && append!(written_inputs, writetojob(f, job, elkinputs; kwargs...))
        # i = length(abiinputs) + 1
        # while i <= length(inputs(job))
        #     i += writetojob(f, job, inputs(job)[i])
        # end
        for i in inputs(job)
            if i ∉ written_inputs
                append!(written_inputs, writetojob(f, job, i; kwargs...))
            end
        end
    end
end

function read_job_line(line)
    line = strip(line)
    line = replace(line, ['>', '<'] => " ")
    if line[1] == '#'
        run = false
        line = strip(line[2:end], '#')
    else
        run = true
    end
    spl = strip_split(line)

    input  = spl[end-1]
    output = spl[end]
    spl    = spl[1:end-2]
    exec_and_flags = Pair{String, Vector{SubString}}[]
    #TODO This is not really nice, we don't handle execs that are unparseable...
    #     Not sure how we can actually do this
    for s in spl
        if any(occursin.(allexecs(), (s,)))
            push!(exec_and_flags, s => SubString[])
        elseif !isempty(exec_and_flags)
        # else
            push!(last(exec_and_flags[end]), s)
        end
    end

    execs = Exec[]
    for (e, flags) in exec_and_flags
        dir, efile = splitdir(e)
        if occursin("mpi", e)
            push!(execs, Exec(efile, dir, parse_mpi_flags(flags)))
        elseif efile == "wannier90.x"
            push!(execs, Exec(efile, dir, parse_wan_execflags(flags)))
        elseif any(occursin.(QE_EXECS, (efile,))) && !occursin("pw2wannier90", efile)
            push!(execs, Exec(efile, dir, parse_qeexecflags(flags)))
        elseif any(occursin.(ELK_EXECS, (efile,)))
	        input = "elk.in"
	        output = "elk.out"
            push!(execs, Exec(efile, dir))
        end
    end
    return execs, input, output, run
end

# TODO: make this work again
# function read_job_filenames(job_file::String)
#     input_files = String[]
#     output_files = String[]
#     open(job_file, "r") do f
#         readline(f)
#         while !eof(f)
#             line = readline(f)
#             if isempty(line)
#                 continue
#             end
#             if occursin(".x", line)
#                 runcommand, exec, input, output, run = read_job_line(line)
#                 !in(input,  input_files)  && push!(input_files,  input)
#                 !in(output, output_files) && push!(output_files, output)
#             end
#         end
#     end
#     return input_files, output_files
# end

function read_job_inputs(job_file::String)
    dir = splitdir(job_file)[1]
    name   = ""
    header = Vector{String}()
    inputs     = DFInput[]
    structures = AbstractStructure[]
    open(job_file, "r") do f
        readline(f)
        while !eof(f)
            line = readline(f)
            if line == ""
                continue
            end
            if has_parseable_exec(line)
                execs, inputfile, output, run = read_job_line(line)
                inpath = joinpath(dir, inputfile)
                if !ispath(inpath)
                    input = (nothing, nothing)
                else
                    calccommand = getfirst(isparseable, execs)
                    input = calccommand != nothing ? inputparser(calccommand)(inpath, execs=execs, run=run) : (nothing, nothing)
                end
                if input != (nothing, nothing)
                    id = findall(x -> infilename(x) == inputfile, inputs)
                    if !isempty(id) #this can only happen for stuff that needs to get preprocessed
                        merge!(flags(input[1]), flags(inputs[id[1]]))
                        inputs[id[1]] = input[1]
                        # structures[id[1]] = input[2]
                    else
	                    if isa(input[1], Vector)
	                        append!(inputs, input[1])
                        else
	                        push!(inputs, input[1])
                        end
                        if input[2] != nothing
                            push!(structures, input[2])
                        end
                    end
                end
	        elseif occursin("#SBATCH", line)
	                if occursin("-J", line)
	                    name = split(line)[end]
	                else
	                    push!(header, line)
	                end
	        else
	            push!(header, line)
	        end
        end
    end
    if isempty(structures)
	    error("Something went wrong and no valid structures could be read from input files.")
    end
    outstruct = mergestructures(structures)
    return name, header, inputs, outstruct
end

#---------------------------END GENERAL SECTION-------------------#

function expr2file(filename::String, expression::Expr)
    eq        = Symbol("=")
    lines     = readlines(filename)
    new_lines = String[]
    found     = false

    if expression.head != eq
        error("For now only writing of assignment expressions is possible.")
    end

    lhs = expression.args[1]
    rhs = expression.args[2]

    for line in lines
        if line == ""
            continue
        end

        expr = Meta.parse(line)
        if typeof(expr) == Nothing
            continue
        end
        if expr.head != eq
            continue
        end

        lhs_t = expr.args[1]
        rhs_t = expr.args[2]

        if lhs_t == lhs
            found = true
            push!(new_lines, "$(:($lhs = $rhs))")
        else
            push!(new_lines, "$expr")
        end
    end

    open(filename, "w") do f
        for line in new_lines
            write(f, line * "\n")
        end
        if !found
            write(f, "$expression\n")
        end
    end
end

function rm_expr_lhs(filename, lhs)
    lines       = readlines(filename)
    write_lines = String[]
    ind_2_rm    = 0

    for line in lines
        lhs_t = Meta.parse(line).args[1]
        if lhs_t == lhs
            continue
        else
            push!(write_lines, line)
        end
    end

    open(filename, "w") do f
        for line in write_lines
            write(f,line * "\n")
        end
    end
end

function write_cell(f::IO, cell::AbstractMatrix)
    @assert size(cell) == (3, 3) "writing cell only allows 3x3 matrices!"
    write(f, "$(cell[1, 1]) $(cell[1, 2]) $(cell[1, 3])\n")
    write(f, "$(cell[2, 1]) $(cell[2, 2]) $(cell[2, 3])\n")
    write(f, "$(cell[3, 1]) $(cell[3, 2]) $(cell[3, 3])\n")
end

"LOL this absolutely is impossible to do for QE"
function writeabortfile(job::DFJob, input::DFInput{QE})
    abortpath = joinpath(job.local_dir,"$(flag(input, :prefix)[2:end-1]).EXIT")
    open(abortpath, "w") do f
        write(f, " \n")
    end
    return abortpath
end

function read_cutoffs_from_pseudofile(file::AbstractString)
    ecutwfc = 0.0
    ecutrho = 0.0
    open(file, "r") do f
        line = readline(f)
        i = 1
        while i < 100 #semi arbitrary cutoff to amount of lines read
            line = readline(f)
            if occursin("Suggested minimum cutoff for wavefunctions:", line)
                ecutwfc = parse(Float64, split(line)[end-1])
                ecutrho = parse(Float64, split(readline(f))[end-1])
                break
            end
            i += 1
        end
    end
    return ecutwfc, ecutrho
end
