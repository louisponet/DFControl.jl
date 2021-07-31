function parse_file(filename::AbstractString, parse_funcs::Vector{<:Pair{String}};
                    extra_parse_funcs::Vector{<:Pair} = Pair{String,Function}[])
    out = Dict{Symbol,Any}()
    open(filename, "r") do f
        while !eof(f)
            line = strip(readline(f))
            if isempty(line)
                continue
            end
            for pf in (parse_funcs, extra_parse_funcs)
                func = getfirst(x -> occursin(x[1], line), pf)
                func === nothing && continue
                try
                    func[2](out, line, f)
                catch 
                    @warn "File corruption or parsing error detected executing parse function $(func[2]).\nTrying to continue smoothly."
                end
            end
        end
    end
    return out
end

include("qe/fileio.jl")
include("abinit/fileio.jl")
include("wannier90/fileio.jl")
include("elk/fileio.jl")
include("vasp/fileio.jl")

#--------------------Used by other file processing------------------#
function parse_k_line(line, T)
    line = replace(replace(line, ")" => " "), "(" => " ")
    splt = split(line)
    k1   = parse(T, splt[4])
    k2   = parse(T, splt[5])
    k3   = parse(T, splt[6])
    w    = parse(T, splt[end])
    return (v = Vec3([k1, k2, k3]), w = w)
end

function write_flag_line(f, flag, data, seperator = "=", i = "")
    flagstr = string(flag)
    if flagstr[end-1] == '_' && tryparse(Int, string(flagstr[end])) != nothing
        flagstr = flagstr[1:end-2] * "($(flagstr[end]))"
    end
    write(f, "  $flagstr$i $seperator ")

    if typeof(data) <: Array
        if length(data) % 3 == 0 && eltype(data) != Int
            write(f, "\n")
            for i in 1:3:length(data)
                write(f, "  $(data[i]) $(data[i + 1]) $(data[i + 2])\n")
            end
        else
            write(f, "  $(data[1])")
            for x in data[2:end]
                write(f, " $x")
            end
            write(f, "\n")
        end

    else #this should work for anything singular valued data such as bools, ''s and other types
        write(f, "$data\n")
    end
end

function parse_flag_val(val, T = Float64)
    if T == String
        return val
    end
    if occursin("d", val) && T <: Number
        val = replace(val, "d" => "e")
    end

    val = replace(val, "." => "")
    t = parse.(eltype(T), split(lowercase(val)))
    #deal with abinit constants -> all flags that are read which are not part of the abi[:structure] get cast into the correct atomic units!
    if length(t) > 1 && typeof(t[end]) == Symbol
        t = t[1:end-1] .* abi_conversions[t[end]]
    end

    return length(t) == 1 ? t[1] : typeof(t) <: Vector{Real} ? convert.(T, t) : t
end

function write_data(f, data)
    if typeof(data) <: Vector{Vector{Float64}} || typeof(data) <: Vector{NTuple{4,Float64}} #k_points
        for x in data
            for y in x
                write(f, " $y")
            end
            write(f, "\n")
        end
    elseif typeof(data) <: Vector{Int} || typeof(data) <: NTuple{6,Int}
        for x in data
            write(f, " $x")
        end
        write(f, "\n")
    elseif typeof(data) <: Matrix
        im, jm = size(data)
        for i in 1:im
            for j in 1:jm
                write(f, " $(data[i, j])")
            end
            write(f, "\n")
        end
    end
end

Base.length(::Type{<:Number}) = 1

function Base.parse(::Type{NamedTuple{names,types}},
                    spl::Vector{<:AbstractString}) where {names,types}
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
Base.parse(::Type{T}, s::AbstractString) where {T<:NamedTuple} = parse(T, split(s))

function Base.parse(::Type{Point{N,T}}, spl::Vector{<:AbstractString}) where {N,T}
    @assert N == length(spl)
    return Point{N,T}(parse.(T, spl))
end
Base.parse(::Type{T}, s::AbstractString) where {T<:Point} = parse(T, split(s))

#---------------------------BEGINNING GENERAL SECTION-------------------#
#Incomplete: only works with SBATCH right now
function write_job_name(f, job::DFJob)
    return write(f, "#SBATCH -J $(job.name) \n")
end

function write_job_header(f, job::DFJob)
    job_header = isempty(job.header) ? getdefault_jobheader() : job.header
    for line in job_header
        if occursin("\n", line)
            write(f, line)
        else
            write(f, line * "\n")
        end
    end
end

function writetojob(f, job, calculations::Vector{DFCalculation{Abinit}}; kwargs...)
    abinit_jobfiles = write_abi_datasets(calculations, job.local_dir; kwargs...)
    abifiles = String[]
    num_abi = 0
    for (filename, pseudos, runcommand) in abinit_jobfiles
        push!(abifiles, filename)
        file, ext = splitext(filename)
        write(f,
              "$runcommand << !EOF\n$filename\n$(file * ".out")\n$(job.name * "_Xi$num_abi")\n$(job.name * "_Xo$num_abi")\n$(job.name * "_Xx$num_abi")\n")
        for pp in pseudos
            write(f, "$pp\n")
        end
        write(f, "!EOF\n")
        num_abi += 1
    end
    return abifiles
end

function writetojob(f, job, calculations::Vector{DFCalculation{Elk}}; kwargs...)
    save(calculations, job.structure; kwargs...)
    should_run = any(map(x -> x.run, calculations))
    if !should_run
        write(f, "#")
    end
    writeexec.((f,), execs(calculations[1]))
    write(f, "< $(infilename(calculations[1])) > $(outfilename(calculations[1]))\n")
    return calculations
end

writeexec(f, exec::Exec) = write(f, string(exec) * " ")

function writetojob(f, job, calculation::DFCalculation; kwargs...)
    filename   = infilename(calculation)
    should_run = calculation.run
    save(calculation, job.structure; kwargs...)
    if !should_run
        write(f, "#")
    end
    writeexec.((f,), execs(calculation))
    write(f, "< $filename > $(outfilename(calculation))\n")
    return (calculation,)
end

function writetojob(f, job, _calculation::DFCalculation{Wannier90}; kwargs...)
    filename   = infilename(_calculation)
    should_run = _calculation.run
    id         = findfirst(isequal(_calculation), job.calculations)
    seedname   = name(_calculation)

    nscf_calc = getfirst(x -> isnscf(x), job.calculations)
    if nscf_calc !== nothing
        runexec = nscf_calc.execs
        # For elk the setup necessary for the wan_calc needs to be done before writing the wan calculation
        # because it's inside elk.in
        if package(nscf_calc) == QE
            pw2wancalculation = qe_generate_pw2wancalculation(_calculation, nscf_calc,
                                                              runexec)
            preprocess = pop!(flags(_calculation), :preprocess)
            wannier_plot = pop!(flags(_calculation), :wannier_plot, nothing)

            if !preprocess || !should_run
                write(f, "#")
            end
            writeexec.((f,), execs(_calculation))
            write(f, "-pp $filename > $(outfilename(_calculation))\n")

            save(_calculation, job.structure; kwargs...)
            writetojob(f, job, pw2wancalculation; kwargs...)
            flags(_calculation)[:preprocess] = preprocess
            wannier_plot !== nothing && (flags(_calculation)[:wannier_plot] = wannier_plot)
        elseif package(nscf_calc) == Elk
            pw2wancalculation = job["elk2wannier"]
        end
    end

    if !should_run
        write(f, "#")
    end
    writeexec.((f,), execs(_calculation))
    write(f, "$filename > $(outfilename(_calculation))\n")
    return (_calculation,)
end

function write_job_preamble(f, job::DFJob)
    if !isempty(job.server_dir)
        dir = splitpath(job.local_dir)[end]
        write(f, "cp -r $(job.local_dir) $(job.server_dir) \n")
        write(f, "cd $(joinpath(job.server_dir, dir)) \n")
    end
end

function write_job_postamble(f, job::DFJob)
    if !isempty(job.server_dir)
        dir = splitpath(job.local_dir)[end]
        write(f, "cp -r $(joinpath(job.server_dir, dir)) $(splitdir(job.local_dir)[1])\n")
        write(f, "rm -r $(joinpath(job.server_dir, dir))\n")
    end
end

"""
    writejobfiles(job::DFJob; kwargs...)

Writes all the calculation files and job file that are linked to a DFJob.
Kwargs will be passed down to various writetojob functions.
"""
function writejobfiles(job::DFJob; kwargs...)
    # rm.(joinpath.(Ref(job.local_dir), searchdir(job.local_dir, ".in")))
    open(joinpath(job.local_dir, "job.tt"), "w") do f
        write(f, "#!/bin/bash\n")
        write_job_name(f, job)
        write_job_header(f, job)
        write_job_preamble(f, job)
        written_calculations = DFCalculation[]
        abicalculations = Vector{DFCalculation{Abinit}}(filter(x -> package(x) == Abinit,
                                                               calculations(job)))
        !isempty(abicalculations) && writetojob(f, job, abicalculations; kwargs...)
        elkcalculations = Vector{DFCalculation{Elk}}(filter(x -> package(x) == Elk,
                                                            calculations(job)))
        !isempty(elkcalculations) &&
            append!(written_calculations, writetojob(f, job, elkcalculations; kwargs...))

        for i in calculations(job)
            if i.run
                rm_outfiles(i)
            end
            if i ∉ written_calculations
                append!(written_calculations, writetojob(f, job, i; kwargs...))
            end
        end
        return write_job_postamble(f, job)
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

    calculation = spl[end-1]
    output = spl[end]
    spl = spl[1:end-2]
    exec_and_flags = Pair{String,Vector{SubString}}[]
    #TODO This is not really nice, we don't handle execs that are unparseable...
    #     Not sure how we can actually do this
    for s in spl
        d, e = splitdir(s)
        if e ∈ allexecs()
            push!(exec_and_flags, s => SubString[])
        elseif !isempty(exec_and_flags)
            # else
            push!(last(exec_and_flags[end]), s)
        end
    end

    execs = Exec[]
    for (e, flags) in exec_and_flags
        dir, efile = splitdir(e)
        dir = replace(dir, "~" => homedir())
        if occursin("pw2wannier90", efile)
            continue
        end
        if occursin("mpi", e)
            push!(execs, Exec(efile, dir, parse_mpi_flags(flags)))
        elseif efile == "wannier90.x"
            push!(execs, Exec(efile, dir, parse_wan_execflags(flags)))
        elseif any(occursin.(QE_EXECS, (efile,)))
            push!(execs, Exec(efile, dir, parse_qe_execflags(flags)))
        elseif any(occursin.(ELK_EXECS, (efile,)))
            calculation = "elk.in"
            output = "elk.out"
            push!(execs, Exec(exec=efile, dir=dir))
        else
            push!(execs, Exec(efile, dir, parse_generic_flags(flags)))
        end
    end
    return execs, calculation, output, run
end

function parse_generic_flags(flags::Vector{<:SubString})
    out = ExecFlag[]
    f = :none
    c = 1
    for s in flags
        if s[1] == '-'
            c = count(isequal('-'), s)
            f = Symbol(strip(s, '-'))
        else
            v = tryparse(Int, s)
            if v === nothing
                v = tryparse(Float64, s)
                if v === nothing
                    v = s
                end
            end
            push!(out, ExecFlag(f => v, c))
        end
    end
    return out
end

# TODO: make this work again
# function read_job_filenames(job_file::String)
#     calculation_files = String[]
#     output_files = String[]
#     open(job_file, "r") do f
#         readline(f)
#         while !eof(f)
#             line = readline(f)
#             if isempty(line)
#                 continue
#             end
#             if occursin(".x", line)
#                 runcommand, exec, calculation, output, run = read_job_line(line)
#                 !in(calculation,  calculation_files)  && push!(calculation_files,  calculation)
#                 !in(output, output_files) && push!(output_files, output)
#             end
#         end
#     end
#     return calculation_files, output_files
# end

function read_job_calculations(job_file::String)
    dir = splitdir(job_file)[1]
    name = ""
    header = Vector{String}()
    calculations = DFCalculation[]
    structures = AbstractStructure[]
    serverdir = ""
    open(job_file, "r") do f
        readline(f)
        while !eof(f)
            line = readline(f)
            if line == ""
                continue
            end
            if has_parseable_exec(line)
                execs, calculationfile, output, run = read_job_line(line)
                inpath = joinpath(dir, calculationfile)
                if !ispath(inpath)
                    calculation = (nothing, nothing)
                else
                    calccommand = getfirst(isparseable, execs)
                    calculation = calccommand != nothing ?
                                  calculationparser(calccommand)(inpath; execs = execs,
                                                                 run = run) :
                                  (nothing, nothing)
                end
                if calculation[1] !== nothing
                    calculation[1].outfile = output
                    calculation[1].infile = calculationfile
                end
                if calculation != (nothing, nothing)
                    id = findall(x -> infilename(x) == calculationfile, calculations)
                    if !isempty(id) #this can only happen for stuff that needs to get preprocessed
                        merge!(flags(calculation[1]), flags(calculations[id[1]]))
                        calculations[id[1]] = calculation[1]
                        # structures[id[1]] = calculation[2]
                    else
                        if isa(calculation[1], Vector)
                            append!(calculations, calculation[1])
                        else
                            push!(calculations, calculation[1])
                        end
                        if calculation[2] != nothing
                            push!(structures, calculation[2])
                        end
                    end
                end
            elseif occursin("#SBATCH", line)
                if occursin("-J", line)
                    name = split(line)[end]
                    name = name == "-J" ? "noname" : name
                else
                    push!(header, line)
                end
            elseif occursin("cp -r", line) && isempty(serverdir)
                serverdir = split(line)[end]
            else
                push!(header, line)
            end
        end
    end
    if isempty(structures)
        error("Something went wrong and no valid structures could be read from calculation files.")
    end
    outstruct = mergestructures(structures)
    return (name = name, header = header, calculations = calculations,
            structure = outstruct, server_dir = serverdir)
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
            write(f, line * "\n")
        end
    end
end

function write_cell(f::IO, cell::AbstractMatrix)
    @assert size(cell) == (3, 3) "writing cell only allows 3x3 matrices!"
    write(f, "$(cell[1, 1]) $(cell[1, 2]) $(cell[1, 3])\n")
    write(f, "$(cell[2, 1]) $(cell[2, 2]) $(cell[2, 3])\n")
    return write(f, "$(cell[3, 1]) $(cell[3, 2]) $(cell[3, 3])\n")
end

"LOL this absolutely is impossible to do for QE"
function writeabortfile(job::DFJob, calculation::DFCalculation{QE})
    abortpath = joinpath(job.local_dir, TEMP_CALC_DIR, "$(job.name).EXIT")
    open(abortpath, "w") do f
        return write(f, " \n")
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

function write_xsf(filename::AbstractString, structure::AbstractStructure)
    open(filename,"w") do f
        write(f, "CRYSTAL\n")
        c = ustrip.(cell(structure)')
        write(f, "PRIMVEC\n")
        write(f, "$(c[1,1]) $(c[1,2]) $(c[1,3])\n")
        write(f, "$(c[2,1]) $(c[2,2]) $(c[2,3])\n")
        write(f, "$(c[3,1]) $(c[3,2]) $(c[3,3])\n")
        write(f, "CONVVEC\n")
        write(f, "$(c[1,1]) $(c[1,2]) $(c[1,3])\n")
        write(f, "$(c[2,1]) $(c[2,2]) $(c[2,3])\n")
        write(f, "$(c[3,1]) $(c[3,2]) $(c[3,3])\n")
        write(f, "PRIMCOORD\n")
        write(f, "$(length(atoms(structure))) 1\n")
        for at in atoms(structure)
            n = element(at).symbol
            p = ustrip.(position_cart(at))
            write(f, "$n $(p[1]) $(p[2]) $(p[3])\n")
        end
    end
end

function writelines(file, lines)
    open(file, "w") do f
        for l in lines
            write(f, l)
            write(f, "\n")
        end
    end
end
