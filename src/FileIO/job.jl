function write_job_header(f, job::Job, environment)
    environment !== nothing && write(f, environment)
    job_header = job.header
    for line in job_header
        if occursin("\n", line)
            write(f, line)
        else
            write(f, line * "\n")
        end
    end
    modules = AbstractString[] 
    for e in map(x->x.exec, job.calculations)
        for m in e.modules
            if !(m ∈ modules)
                push!(modules, m)
            end
        end
    end
    if !isempty(modules)
        write(f, """module load $(join(modules, " "))\n""")
    end
end

# function writetojob(f, job, calculations::Vector{Calculation{Abinit}}, environment; kwargs...)
#     abinit_jobfiles = write_abi_datasets(calculations, job.dir; kwargs...)
#     abifiles = String[]
#     num_abi = 0
#     for (filename, pseudos, runcommand) in abinit_jobfiles
#         push!(abifiles, filename)
#         file, ext = splitext(filename)
#         write(f,
#               "$runcommand << !EOF\n$filename\n$(file * ".out")\n$(job.name * "_Xi$num_abi")\n$(job.name * "_Xo$num_abi")\n$(job.name * "_Xx$num_abi")\n")
#         for pp in pseudos
#             write(f, "$pp\n")
#         end
#         write(f, "!EOF\n")
#         num_abi += 1
#     end
#     return abifiles
# end

function writetojob(f, job, calculations::Vector{Calculation{Elk}}, environment; kwargs...)
    save(calculations, job.structure, job.dir; kwargs...)
    should_run = any(map(x -> x.run, calculations))
    if !should_run
        write(f, "#")
    end
    writeexec(f, calculations[1].exec, environment)
    write(f, "< $(calculations[1].infile) > $(calculations[1].outfile)\n")
    return calculations
end

function writeexec(f, exec::Exec, environment)
    if exec.parallel
        @assert environment !== nothing "Exec $(exec.exec) flagged as parallel but no valid environment was supplied."
        write(f, environment.MPI_command * " ")
    end
    write(f, string(exec) * " ")
end

function writetojob(f, job, calculation::Calculation, environment; kwargs...)
    filename   = calculation.infile
    should_run = calculation.run
    save(calculation, job.structure, joinpath(job, calculation.infile); kwargs...)
    if !should_run
        write(f, "#")
    end
    writeexec(f, calculation.exec, environment)
    write(f, "< $filename > $(calculation.outfile)\n")
    return (calculation,)
end

function writetojob(f, job, _calculation::Calculation{Wannier90}, environment; kwargs...)
    filename   = _calculation.infile
    should_run = _calculation.run
    id         = findfirst(isequal(_calculation), job.calculations)
    seedname   = _calculation.name

    nscf_calc = getfirst(x -> Calculations.isnscf(x), job.calculations)
    if nscf_calc !== nothing
        runexec = nscf_calc.exec
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
            writeexec(f, _calculation.exec, environment)
            write(f,
                  "-pp $filename > $(_calculation.outfile)\n")

            save(_calculation, job.structure, joinpath(job, _calculation.infile); kwargs...)
            writetojob(f, job, pw2wancalculation, environment; kwargs...)
            _calculation[:preprocess] = preprocess
            wannier_plot !== nothing && (_calculation[:wannier_plot] = wannier_plot)
        elseif eltype(nscf_calc) == Elk
            pw2wancalculation = job["elk2wannier"]
        end
    end

    if !should_run
        write(f, "#")
    end
    writeexec(f, _calculation.exec, environment)
    write(f, "$filename > $(_calculation.outfile)\n")
    return (_calculation,)
end

#TODO: This should take scratch dir of server

function write_job_preamble(f, job::Job) end
function write_job_postamble(f, job::Job) end
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
    write_job_files(job::Job; kwargs...)

Writes all the calculation files and job file that are linked to a Job.
Kwargs will be passed down to various writetojob functions.
"""
function write_job_files(job::Job, environment; kwargs...)

    if any(x ->eltype(x) == QE, job.calculations)
        for a in unique(job.structure.atoms)
            write(joinpath(job, "$(a.element.symbol).UPF"), a.pseudo)
        end
    end
            
    open(joinpath(job, "job.tt"), "w") do f
        write(f, "#!/bin/bash\n")
        write(f, "#SBATCH -J $(job.name) \n")
        write_job_header(f, job, environment)
        write_job_preamble(f, job)
        written_calculations = Calculation[]
        abicalculations = filter(x -> eltype(x) == Abinit, job.calculations)
        !isempty(abicalculations) && writetojob(f, job, abicalculations, environment; kwargs...)
        elkcalculations = filter(x -> eltype(x) == Elk, job.calculations)
        !isempty(elkcalculations) &&
            append!(written_calculations, writetojob(f, job, elkcalculations, environment; kwargs...))

        for i in job.calculations
            if i ∉ written_calculations
                append!(written_calculations, writetojob(f, job, i, environment; kwargs...))
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
            push!(execs, Exec(exec=efile, dir=dir, flags=Calculations.parse_mpi_flags(flags)))
        elseif efile == "wannier90.x"
            push!(execs, Exec(exec=efile, dir=dir, flags=Calculations.parse_wan_execflags(flags)))
        elseif any(occursin.(Calculations.QE_EXECS, (efile,)))
            push!(execs, Exec(exec=efile, dir=dir, flags=Calculations.parse_qe_execflags(flags)))
        elseif any(occursin.(Calculations.ELK_EXECS, (efile,)))
            calculation = "elk.in"
            output = "elk.out"
            push!(execs, Exec(; exec = efile, dir = dir))
        else
            push!(execs, Exec(exec=efile, dir=dir, flags=Calculations.parse_generic_flags(flags)))
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

function read_job_script(job_file::String)
    name = ""
    dir = splitdir(job_file)[1]
    header = Vector{String}()
    scratch_dir = ""
    calcs = NamedTuple{(:exec, :infile, :outfile, :run), Tuple{Exec, String, String, Bool}}[]
    open(job_file, "r") do f
        readline(f)
        while !eof(f)
            line = readline(f)
            
            isempty(line) && continue
            
            if occursin("#SBATCH", line)
                if occursin("-J", line) || occursin("--job-name", line)
                    name = split(split(line, "=")[end])[end]
                else
                    push!(header, line)
                end
            elseif occursin("cp -r", line) && isempty(scratch_dir)
                scratch_dir = split(line)[end]
            elseif Calculations.has_parseable_exec(line)
                execs, infile, outfile, run = read_job_line(line)
                inpath = joinpath(dir, infile)
                if ispath(inpath)
                    exec = getfirst(Calculations.isparseable, execs)
                    if exec !== nothing
                        if length(execs) == 1
                            exec.parallel = false
                        end
                        push!(calcs, (exec = exec, infile = joinpath(dir, infile), outfile = joinpath(dir, outfile), run = run))
                    end
                end
            else
                push!(header, line)
            end
        end
    end

    module_line_ids = findall(x -> occursin("module", x), header)
    
    if module_line_ids !== nothing
        module_lines = header[module_line_ids]
        deleteat!(header, module_line_ids)
        modules = map(x -> split(x)[end], module_lines)
    end
    #TODO cleanup
    known_es = Calculations.load_execs()
    execs = filter(x -> !any(y -> y.dir == x.dir && y.exec == x.exec, values(known_es)), unique(map(x->x.exec, calcs)))
    for e in execs
        i = 1
        runnable = Calculations.isrunnable(e)
        while length(e.modules) < length(modules) && !runnable
            push!(e.modules, modules[i])
            runnable = Calculations.isrunnable(e)
            i+=1
        end
        if runnable
            e.name = """$(e.exec)@$(join(splitpath(e.dir)[2:end], "_"))"""
            Calculations.maybe_register(e)
        end
    end
    for c in calcs
        e1 = c.exec
        for (n, e) in merge(Dict([e.name => e for e in execs]), known_es)
            found = false
            if e.dir == e1.dir && e.exec == e1.exec
                e1.modules = e.modules
            end
        end
    end

    runtime = Jobs.environment_from_jobscript(job_file)
    ename = Jobs.environment_name(runtime)
    if ename !== nothing
        environment = ename
        #TODO only works with slurm!!!
        deleteat!(header, findall(x -> x[1:7] == "#SBATCH", header))
        deleteat!(header, findall(x -> x[1:6] == "export", header))
    else
        environment = ""
    end
    return (name, calcs, header, environment) 
end

function parse_calculations(calcs)
    structures = Structure[]
    outcalcs = Calculation[]
    for calc in calcs
        c = calculationparser(calc.exec)(calc.infile; exec = calc.exec, infile = splitpath(calc.infile)[end], outfile = splitpath(calc.outfile)[end], run = calc.run)
        if c[2] !== nothing
            push!(structures, c[2])
        end
        if c[1] !== nothing
            push!(outcalcs, c[1])
        end
    end
    if !isempty(structures)
        structure = Structures.mergestructures(structures)
    else
        structure = Structure()
        @warn "No valid structures could be read from calculation files."
    end
    return (calculations = outcalcs, structure=structure)
end
