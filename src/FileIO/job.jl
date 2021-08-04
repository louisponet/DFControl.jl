function write_job_header(f, job::Job)
    job_header = job.header
    for line in job_header
        if occursin("\n", line)
            write(f, line)
        else
            write(f, line * "\n")
        end
    end
end

function writetojob(f, job, calculations::Vector{Calculation{Abinit}}; kwargs...)
    abinit_jobfiles = write_abi_datasets(calculations, job.dir; kwargs...)
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

function writetojob(f, job, calculations::Vector{Calculation{Elk}}; kwargs...)
    save(calculations, job.structure; kwargs...)
    should_run = any(map(x -> x.run, calculations))
    if !should_run
        write(f, "#")
    end
    writeexec.((f,), calculations[1].execs)
    write(f, "< $(calculations[1].infile) > $(calculations[1].outfile)\n")
    return calculations
end

writeexec(f, exec::Exec) = write(f, string(exec) * " ")

function writetojob(f, job, calculation::Calculation; kwargs...)
    filename   = calculation.infile
    should_run = calculation.run
    save(calculation, job.structure; kwargs...)
    if !should_run
        write(f, "#")
    end
    writeexec.((f,), calculation.execs)
    write(f, "< $filename > $(calculation.outfile)\n")
    return (calculation,)
end

function writetojob(f, job, _calculation::Calculation{Wannier90}; kwargs...)
    filename   = _calculation.outfile
    should_run = _calculation.run
    id         = findfirst(isequal(_calculation), job.calculations)
    seedname   = _calculation.name

    nscf_calc = getfirst(x -> Calculations.isnscf(x), job.calculations)
    if nscf_calc !== nothing
        runexec = nscf_calc.execs
        # For elk the setup necessary for the wan_calc needs to be done before writing the wan calculation
        # because it's inside elk.in
        if eltype(nscf_calc) == QE
            pw2wancalculation = qe_generate_pw2wancalculation(_calculation, nscf_calc,
                                                              runexec)
            preprocess = pop!(_calculation, :preprocess, false)
            wannier_plot = pop!(_calculation, :wannier_plot, nothing)

            if !preprocess || !should_run
                write(f, "#")
            end
            writeexec.((f,), _calculation.execs)
            write(f, "-pp $filename > $(joinpath(_calculation.dir, _calculation.outfile))\n")

            save(_calculation, job.structure; kwargs...)
            writetojob(f, job, pw2wancalculation; kwargs...)
            _calculation[:preprocess] = preprocess
            wannier_plot !== nothing && (_calculation[:wannier_plot] = wannier_plot)
        elseif eltype(nscf_calc) == Elk
            pw2wancalculation = job["elk2wannier"]
        end
    end

    if !should_run
        write(f, "#")
    end
    writeexec.((f,), _calculation.execs)
    write(f, "$filename > $(joinpath(_calculation.dir, _calculation.outfile))\n")
    return (_calculation,)
end

#TODO: This should take scratch dir of server

function write_job_preamble(f, job::Job)
end
function write_job_postamble(f, job::Job)
end
# function write_job_preamble(f, job::Job)
#     if !isempty(job.server_dir)
#         dir = splitpath(job.dir)[end]
#         write(f, "cp -r $(job.dir) $(job.server_dir) \n")
#         write(f, "cd $(joinpath(job.server_dir, dir)) \n")
#     end
# end

# function write_job_postamble(f, job::Job)
#     if !isempty(job.server_dir)
#         dir = splitpath(job.dir)[end]
#         write(f, "cp -r $(joinpath(job.server_dir, dir)) $(splitdir(job.dir)[1])\n")
#         write(f, "rm -r $(joinpath(job.server_dir, dir))\n")
#     end
# end

"""
    writejobfiles(job::Job; kwargs...)

Writes all the calculation files and job file that are linked to a Job.
Kwargs will be passed down to various writetojob functions.
"""
function writejobfiles(job::Job; kwargs...)
    # rm.(joinpath.(Ref(job.dir), searchdir(job.dir, ".in")))
    open(joinpath(job.dir, "job.tt"), "w") do f
        write(f, "#!/bin/bash\n")
        write(f, "#SBATCH -J $(job.name) \n")
        write_job_header(f, job)
        write_job_preamble(f, job)
        written_calculations = Calculation[]
        abicalculations = filter(x -> eltype(x) == Abinit, job.calculations)
        !isempty(abicalculations) && writetojob(f, job, abicalculations; kwargs...)
        elkcalculations = filter(x -> eltype(x) == Elk, job.calculations)
        !isempty(elkcalculations) &&
            append!(written_calculations, writetojob(f, job, elkcalculations; kwargs...))

        for i in job.calculations
            if i.run
                rm.(Calculations.outfiles(i))
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
        if e ∈ Calculations.allexecs()
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
            push!(execs, Exec(efile, dir, Calculations.parse_mpi_flags(flags)))
        elseif efile == "wannier90.x"
            push!(execs, Exec(efile, dir, Calculations.parse_wan_execflags(flags)))
        elseif any(occursin.(Calculations.QE_EXECS, (efile,)))
            push!(execs, Exec(efile, dir, Calculations.parse_qe_execflags(flags)))
        elseif any(occursin.(Calculations.ELK_EXECS, (efile,)))
            calculation = "elk.in"
            output = "elk.out"
            push!(execs, Exec(exec=efile, dir=dir))
        else
            push!(execs, Exec(efile, dir, Calculations.parse_generic_flags(flags)))
        end
    end
    return execs, calculation, output, run
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

function calculationparser(exec::Exec)
    if Calculations.is_qe_exec(exec)
        qe_read_calculation
    elseif Calculations.is_wannier_exec(exec)
        wan_read_calculation
    elseif Calculations.is_elk_exec(exec)
        elk_read_calculation
    end
end


function read_job_calculations(job_file::String)
    dir = splitdir(job_file)[1]
    name = ""
    header = Vector{String}()
    cs = Calculation[]
    structures = Structure[]
    scratch_dir = ""
    open(job_file, "r") do f
        readline(f)
        while !eof(f)
            line = readline(f)
            if line == ""
                continue
            end
            if Calculations.has_parseable_exec(line)
                execs, cfile, output, run = read_job_line(line)
                inpath = joinpath(dir, cfile)
                if !ispath(inpath)
                    c = (nothing, nothing)
                else
                    command = getfirst(Calculations.isparseable, execs)
                    c = command != nothing ?
                                  calculationparser(command)(inpath; execs = execs,
                                                                 run = run) :
                                  (nothing, nothing)
                end
                if c[1] !== nothing
                    c[1].outfile = output
                    c[1].infile = cfile
                end
                if c != (nothing, nothing)
                    id = findall(x -> x.infile == cfile, cs)
                    if !isempty(id) #this can only happen for stuff that needs to get preprocessed
                        merge!(c[1].flags, cs[id[1]].flags)
                        cs[id[1]] = c[1]
                    else
                        if isa(c[1], Vector)
                            append!(cs, c[1])
                        else
                            push!(cs, c[1])
                        end
                        if c[2] != nothing
                            push!(structures, c[2])
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
            elseif occursin("cp -r", line) && isempty(scratch_dir)
                scratch_dir = split(line)[end]
            else
                push!(header, line)
            end
        end
    end
    if !isempty(structures)
        structure = Structures.mergestructures(structures)
    else
        structure = Structure()
        @warn "No valid structures could be read from calculation files."
    end
    return (;name, header, calculations=cs, structure)
end

