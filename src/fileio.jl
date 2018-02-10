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

    length(t) == 1 ? t[1] : typeof(t) <: Array{Real,1} ? convert.(T,t) : t
end

#---------------------------BEGINNING GENERAL SECTION-------------------#
"""
    write_input(df_input::DFInput, filename::String=df_input.filename)

Writes the input file for a DFInput.Backend of DFInput decides what writing function is called.
"""
function write_input(df_input::DFInput, filename::String=df_input.filename)
    if typeof(df_input) == QEInput
        write_qe_input(df_input, filename)
    elseif typeof(df_input) == WannierInput
        write_wannier_input(df_input, filename)
    elseif typeof(df_input) == AbinitInput
        write_abi_input(df_input, filename)
    end
end

#Incomplete: only works with SBATCH right now
function write_job_name(job::DFJob, f)
    write(f, "#SBATCH -J $(job.name) \n")
end

function write_job_header(job::DFJob, f)
    job_header = job.header == "" ? get_default_job_header() : job.header
    for line in job_header
        if contains(line, "\n")
            write(f, line)
        else
            write(f, line * "\n")
        end
    end
end

"""
    write_job_files(job::DFJob)

Writes all the input files and job file that are linked to a DFJob.
"""
function write_job_files(job::DFJob)
    files_to_remove = search_dir(job.local_dir, ".in")
    new_filenames   = String[]
    num_abi         = 0
    open(job.local_dir * "job.tt", "w") do f
        write(f, "#!/bin/bash\n")
        write_job_name(job, f)
        write_job_header(job, f)
        i = 1
        while i <= length(job.calculations)
            calculation = job.calculations[i]
            run_command = calculation.run_command
            filename    = calculation.filename
            should_run  = calculation.run
            if typeof(calculation) == WannierInput
                write_input(calculation, job.local_dir * filename)
                if calculation.preprocess
                    run_command *= " -pp"
                end
                if !should_run
                    write(f, "#$run_command $(filename[1:end-4])\n")
                else
                    write(f, "$run_command $(filename[1:end-4])\n")
                end
                i += 1
            elseif typeof(calculation) == AbinitInput
                abinit_inputs     = Array{AbinitInput,1}(filter(x -> typeof(x) == AbinitInput, job.calculations))
                i += length(abinit_inputs)
                abinit_jobfiles   = write_abi_datasets(abinit_inputs, job.local_dir)
                for (filename, pseudos, run_command) in abinit_jobfiles
                    file, ext = splitext(filename)
                    write(f, "$run_command << !EOF\n$filename\n$(file * ".out")\n$(job.name * "_Xi$num_abi")\n$(job.name * "_Xo$num_abi")\n$(job.name * "_Xx$num_abi")\n")
                    for pp in pseudos
                        write(f, "$pp\n")
                    end
                    write(f, "!EOF\n")
                    num_abi += 1
                end

            elseif typeof(calculation) == QEInput
                exec = calculation.exec
                write_input(calculation, job.local_dir * filename)
                if !should_run
                    write(f, "#$run_command $exec <$filename> $(split(filename,".")[1]).out \n")
                else
                    write(f, "$run_command $exec <$filename> $(split(filename,".")[1]).out \n")
                end
                i += 1
            end
            push!(new_filenames, filename)
        end
    end

    for file in files_to_remove
        if !in(file, new_filenames)
            rm(job.local_dir * file)
        end
    end
    return new_filenames
end


function read_command_line(line)
    if typeof(line) <: String
        line = split(line)
    end
    i = 0
    for (j, s) in enumerate(line)
        if contains(s, ".x")
            i = j
            break
        elseif contains(s, "abinit")
            i = j
            break
        end
    end
    run_command = prod([s * " " for s in line[1:i]])
#TODO think about a better way of handling wannier90, probably with a fixed set of rules. Only 1 wannier90 input, which preprocesses or not..
    if contains(line[i + 1], "-pp") #this is for wannier90
        run_command *= line[i + 1]
        i += 1
    end
    return i, run_command
end


"""
    read_job_file(job_file::String)

Reads and returns the name, input files, run_commands and whether or not they need to be commented out.
All files that are read contain "in".
This reads QE and wannier90 inputs for now.
"""
function read_job_file(job_file::String)
    data = Dict{Symbol,Any}()
    data[:name]         = ""
    data[:header]       = Array{String,1}()
    data[:input_files]  = Array{String,1}()
    data[:output_files] = Array{String,1}()
    data[:run_commands] = Array{String,1}()
    data[:should_run]   = Array{Bool,1}()
    open(job_file, "r") do f
        readline(f)
        while !eof(f)
            line = readline(f)
            if line == ""
                continue
            end
            if contains(line, ".x ")
                if contains(line, "#")
                    push!(data[:should_run], false)
                    line = line[2:end]
                else
                    push!(data[:should_run], true)
                end

                s_line        = split(line)
                i, run_command = read_command_line(s_line)
                push!(data[:run_commands], run_command)
                #handles QE and Wannier.
                in_out = strip_split(prod(s_line[i + 1:end]), ">")
                if length(in_out) == 2
                    push!(data[:input_files],  strip(in_out[1], '<'))
                    push!(data[:output_files], in_out[2])
                else
                    input_file = strip(in_out[1], '<')
                    push!(data[:input_files], input_file)
                    if contains(run_command, "wannier90.x")
                        push!(data[:output_files], input_file * ".wout")
                    end
                end
            elseif contains(line, "abinit ")
                data[:abinit_pseudos] = Array{String,1}()
                s_line         = split(line)
                i, run_command = read_command_line(s_line)
                push!(data[:run_commands], strip(run_command, '#'))
                if contains(line, "!EOF")
                    push!(data[:input_files],  strip(readline(f), '#'))
                    push!(data[:output_files], strip(readline(f), '#'))
                    push!(data[:should_run], true)
                    line = readline(f)
                    while !contains(line, "EOF")
                        if contains(line, ".xml")
                            push!(data[:abinit_pseudos], strip(line, '#'))
                        end
                        line = readline(f)
                    end
                end

                #this is reading the sbatch lines
            elseif contains(line, "#SBATCH")
                if contains(line, "-J")
                    data[:name] = split(line)[end]
                else
                    push!(data[:header], line)
                end
            else
                push!(data[:header], line)
            end
        end
    end
    return data
end

#Incomplete: because QE is stupid we have to find a way to distinguish nscf and bands outputs hardcoded.
# function read_output(filename::string, args...)
#   open(filename,"r") do f
#     while !eof(f)
#       line = readline(f)

#       if contains(line,"self-consistent calculation")
#         return
#       elseif contains(line,"band structure calculation") && !contains(filename,"nscf")
#         return read_qe_bands_file(filename,args...)
#       elseif contains(line, "


#---------------------------END GENERAL SECTION-------------------#

#Incomplete: we should probably handle writing an array of expressions as well!
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
