include("qe/fileio.jl")
include("abinit/fileio.jl")
include("wannier90/fileio.jl")


#--------------------Used by other file processing------------------#
function parse_k_line(line, T)
    splt = split(line)
    k1   = parse(T, splt[5])
    k2   = parse(T, splt[6])
    k3   = parse(T, splt[7][1:1:end-2])
    return [k1, k2, k3]
end

function write_flag_line(f, flag, data, seperator="=", i="")
    write(f,"  $flag$i $seperator ")

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

    if contains(val, "d")
        val = replace(val, "d", "e")
    end

    val = strip(val, '.')
    t = convert.(T, parse.(split(lowercase(val))))
    #deal with abinit constants -> all flags that are read which are not part of the abi[:structure] get cast into the correct atomic units!
    if length(t) > 1 && typeof(t[end]) == Symbol
        t = t[1:end-1] .* abi_conversions[t[end]]
    end

    length(t) == 1 ? t[1] : typeof(t) <: Vector{Real} ? convert.(T,t) : t
end

#---------------------------BEGINNING GENERAL SECTION-------------------#
#Incomplete: only works with SBATCH right now
function write_job_name(f, job::DFJob)
    write(f, "#SBATCH -J $(job.name) \n")
end

function write_job_header(f, job::DFJob)
    job_header = job.header == "" ? get_default_job_header() : job.header
    for line in job_header
        if contains(line, "\n")
            write(f, line)
        else
            write(f, line * "\n")
        end
    end
end

function writetojob(f, job, inputs::Vector{AbinitInput})
    abinit_jobfiles   = write_abi_datasets(inputs, job.local_dir)
    abifiles = String[]
    num_abi = 0
    for (filename, pseudos, run_command) in abinit_jobfiles
        push!(abifiles, filename)
        file, ext = splitext(filename)
        write(f, "$run_command << !EOF\n$filename\n$(file * ".out")\n$(job.name * "_Xi$num_abi")\n$(job.name * "_Xo$num_abi")\n$(job.name * "_Xx$num_abi")\n")
        for pp in pseudos
            write(f, "$pp\n")
        end
        write(f, "!EOF\n")
        num_abi += 1
    end
    return abifiles
end

function writeexec(f, exec::Exec)
    direxec = joinpath(exec.dir, exec.exec)
    write(f, "$direxec")
    for (flag, value) in exec.flags
        write(f," -$flag $value")
    end
    write(f, " ")
end

function writetojob(f, job, input::QEInput)
    exec        = input.exec
    run_command = input.run_command
    filename    = input.filename
    should_run  = input.run
    write_input(input, job.structure, job.local_dir * filename)
    if !should_run
        write(f, "#")
    end
    writeexec(f, run_command)
    writeexec(f, exec)
    write(f, "< $filename > $(split(filename,".")[1]).out\n")
    return 1
end

function writetojob(f, job, input::WannierInput)
    run_command = input.run_command
    filename    = input.filename
    should_run  = input.run
    exec        = input.exec
    # run_command = join([run_command[1], exec[1]], " ")
    #cleanup ugly
    id = findfirst(job.calculations, input)
    seedname = splitext(filename)[1]

    pw2wanid = findfirst(x -> x.control_blocks[1].name == :inputpp && contains(x.exec.exec, "pw2wannier90"), job.calculations[id+1:end])+id
    pw2wan   = job.calculations[pw2wanid]
    stfls!(pw2wan, :seedname => "'$(splitext(input.filename)[1])'")

    if !pw2wan.run
        write(f, "#")
    end
    writeexec(f, run_command)
    writeexec(f, exec)
    write(f, "-pp $(filename[1:end-4]).win > $(filename[1:end-4]).wout\n")

    write_input(input, job.structure, job.local_dir * filename)
    writetojob(f, job, pw2wan)

    if !should_run
        write(f, "#")
    end
    writeexec(f, run_command)
    writeexec(f, exec)
    write(f, "$(filename[1:end-4]).win > $(filename[1:end-4]).wout\n")
    return 2
end
"""
    write_job_files(job::DFJob)

Writes all the input files and job file that are linked to a DFJob.
"""
function write_job_files(job::DFJob)
    files_to_remove = search_dir(job.local_dir, ".in")
    new_filenames   = getfield.(job.calculations, :filename)
    open(job.local_dir * "job.tt", "w") do f
        write(f, "#!/bin/bash\n")
        write_job_name(f, job)
        write_job_header(f, job)
        abiinputs = Vector{AbinitInput}(filter(x -> typeof(x) == AbinitInput, job.calculations))
        !isempty(abiinputs) && push!(new_filenames, writetojob(f, job, abiinputs)...)
        i = length(abiinputs) + 1
        while i <= length(job.calculations)
            i += writetojob(f, job, job.calculations[i])
        end
    end

    for file in files_to_remove
        if !in(file, new_filenames)
            rm(job.local_dir * file)
        end
    end
    return new_filenames
end

function read_job_line(line)
    line = strip(line)
    line = replace(line, ['>', '<'], " ")

    if line[1] == '#'
        run = false
        line = line[2:end]
    else
        run = true
    end

    spl          = strip_split(line)
    run_commands = Exec[]
    t_runcommand = ""
    t_flags      = Dict{Symbol, Any}()
    i = 1
    while i <= length(spl) - 2
        ts = spl[i]
        if ts[1] != '-'
            if t_runcommand == ""
                t_runcommand = ts
            else
                dir, file = splitdir(t_runcommand)
                push!(run_commands, Exec(file, dir, t_flags))
                t_runcommand = ts
                t_flags = Dict{Symbol, Any}()
            end
            i += 1
        else
            if ts == "-pp"
                i += 1
                continue
            end
            flag = Symbol(ts[2:end])
            val  = parse(spl[i + 1])
            t_flags[flag] = val
            i += 2
        end
    end
    dir, file = splitdir(t_runcommand)
    push!(run_commands, Exec(file, dir, t_flags))

    run_command  = length(run_commands) == 1 ? Exec("") : run_commands[1]
    exec         = run_commands[end]
    input        = spl[end-1]
    output       = spl[end]
    return run_command, exec, input, output, run
end
# TODO: make this work again
function read_job_filenames(job_file::String)
    input_files = String[]
    output_files = String[]
    open(job_file, "r") do f
        readline(f)
        while !eof(f)
            line = readline(f)
            if isempty(line)
                continue
            end
            if contains(line, ".x")
                run_command, exec, input, output, run = read_job_line(line)
                !in(input,  input_files)  && push!(input_files,  input)
                !in(output, output_files) && push!(output_files, output)
            end
        end
    end
    return input_files, output_files
end

"""
    read_job_file(job_file::String)

Reads and returns all the relevant information contained in the job input file.
All files that are read contain extension "in".
"""
function read_job_inputs(job_file::String)
    dir = splitdir(job_file)[1]
    name   = ""
    header = Vector{String}()
    inputs     = DFInput[]
    structures = Union{AbstractStructure, Void}[]
    open(job_file, "r") do f
        readline(f)
        while !eof(f)
            line = readline(f)
            if line == ""
                continue
            end
            if contains(line, ".x ")
                run_command, exec, inputfile, output, run = read_job_line(line)
                only_exec = exec.exec
                if only_exec in parseable_qe_execs
                    input, structure = read_qe_input(joinpath(dir, inputfile), run_command=run_command, run=run, exec=exec)
                elseif only_exec == "wannier90.x"
                    input, structure = read_wannier_input(joinpath(dir, splitext(inputfile)[1] * ".win"), run_command=run_command, run=run, exec=exec)
                end
                id = find(x-> x == splitext(inputfile)[1], getindex.(splitext.(getfield.(inputs, :filename)),1))
                if !isempty(id)
                    inputs[id[1]] = input
                    structures[id[1]] = structure
                else
                    push!(inputs, input)
                    push!(structures, structure)
                end
            #Incomplete: Handle abinit in the new way!
            # elseif contains(line, "abinit ")
            #     data[:abinit_pseudos] = Array{String,1}()
            #     s_line         = split(line)
            #     i, run_command = read_job_line(s_line)
            #     push!(data[:run_commands], strip(run_command, '#'))
            #     if contains(line, "!EOF")
            #         push!(data[:input_files],  strip(readline(f), '#'))
            #         push!(data[:output_files], strip(readline(f), '#'))
            #         push!(data[:should_run], true)
            #         line = readline(f)
            #         while !contains(line, "EOF")
            #             if contains(line, ".xml")
            #                 push!(data[:abinit_pseudos], strip(line, '#'))
            #             end
            #             line = readline(f)
            #         end
            #     end
            #
            #     #this is reading the sbatch lines
            elseif contains(line, "#SBATCH")
                if contains(line, "-J")
                    name = split(line)[end]
                else
                    push!(header, line)
                end
            else
                push!(header, line)
            end
        end
    end
    outstruct = merge_structures(structures)
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

        expr = parse(line)
        if typeof(expr) == Void
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
        lhs_t = parse(line).args[1]
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