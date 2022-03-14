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
        # For elk the setup necessary for the wan_calc needs to be done before writing the wan calculation
        # because it's inside elk.in
        if eltype(nscf_calc) == QE
            pw2wancalculation = qe_generate_pw2wancalculation(_calculation, nscf_calc)
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
    write(job::Job, environment::Environment; kwargs...)

Writes all the calculation files and job file that are linked to a Job.
Kwargs will be passed down to various writetojob functions.
"""
function Base.write(job::Job, environment::Environment; kwargs...)

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

    calculation = String(spl[end-1])
    output = String(spl[end])
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
    name::String = ""
    dir = splitdir(job_file)[1]
    header = Vector{String}()
    scratch_dir::String = ""
    calcs = NamedTuple{(:exec, :infile, :outfile, :run), Tuple{Exec, String, String, Bool}}[]
    open(job_file, "r") do f
        readline(f)
        while !eof(f)
            line = readline(f)
            
            isempty(line) && continue
            
            if occursin("#SBATCH", line)
                if occursin("-J", line) || occursin("--job-name", line)
                    name = String(split(split(line, "=")[end])[end])
                else
                    push!(header, line)
                end
            elseif occursin("cp -r", line) && isempty(scratch_dir)
                scratch_dir = String(split(line)[end])
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
    else
        modules = String[]
    end
    #TODO cleanup
    ues = unique(map(x->x.exec, filter(y->y.run, calcs)))
    for e in ues
        t = Database.replacements(e)
        if !isempty(t)
            b = getfirst(x -> e.dir == x.dir && e.exec == x.exec, t)
            if b !== nothing
                e.name = b.name
                e.modules = b.modules
            end
        else
            e.modules = modules
        end
    end
            
    for c in calcs
        e1 = c.exec
        for e in ues
            if e.dir == e1.dir && e.exec == e1.exec
                e1.modules = e.modules
                e1.name = e.name
            end
        end
    end

    runtime = Jobs.environment_from_jobscript(job_file)
    environment = runtime.name
    deleteat!(header, findall(x -> occursin("#SBATCH", x), header))
    deleteat!(header, findall(x -> occursin("export", x), header))
    return (name, calcs, header, environment) 
end

function parse_calculations(calcs)
    structures = Structure[]
    outcalcs = Calculation[]
    for calc in calcs
        exec = Exec(calc[:exec])
        infile = splitpath(calc[:infile])[end]
        if Calculations.is_wannier_exec(exec) && !isempty(outcalcs) && outcalcs[end].infile == infile
            Calculations.set_flags!(outcalcs[end], :preprocess => outcalcs[end].run, print=false)
            empty!(outcalcs[end].exec.flags)
        else
            c = calculationparser(exec)(calc[:contents])
            if c.structure !== nothing
                push!(structures, c.structure)
            end
            push!(outcalcs, Calculation(calc[:name], c.flags, c.data, exec, calc[:run], calc[:infile], calc[:outfile]))
        end
    end
    if !isempty(structures)
        structure = Structures.mergestructures(structures)
    else
        structure = Structure()
        @warn "No valid structures could be read from calculation files."
    end
    Calculations.rm_tmp_flags!.(outcalcs)    
    return (calculations = outcalcs, structure=structure)
end

function Base.write(workflow::Jobs.Workflow, job::Job)
    wdir = mkpath(joinpath(job.dir, ".workflow"))
    qdir = mkpath(joinpath(wdir,"queued"))
    mkpath(joinpath(wdir,"finished"))
    cp(workflow.project_path, joinpath(wdir, "Project.toml"))
    write(joinpath(wdir, "run.jl"),
    """
    cd(@__DIR__)
    using Pkg
    Pkg.instantiate()
    Pkg.develop("DFControl")
    using $(join(workflow.modules, ", "))

    function run_queue(jobdir::String)
        job = DFC.Service.load_job(jobdir)
        logger = DFC.Service.workflow_logger(job)
        DFC.Service.with_logger(logger) do
            qd = DFC.Service.queued_dir(job)
            queued_steps = readdir(qd)
            order = sortperm(parse.(Int, getindex.(splitext.(queued_steps), 1)))
            fd = DFC.Service.finished_dir(job)
            if !ispath(qd) || isempty(queued_steps)
                return
            end
            if !ispath(fd)
                mkpath(fd)
            else
                for f in readdir(fd)
                    rm(joinpath(fd, f))
                end
            end
            for f in queued_steps[order]
                step_file = joinpath(qd, f)
                t = include(step_file)
                results = @eval \$(t)(\$(job))
                DFC.Service.JLD2.jldopen(joinpath(job, ".workflow/results.jld2"), "w+") do j
                    j[splitext(f)[1]] = results
                end
                mv(step_file, joinpath(fd, f); force = true)
                info = state(jobdir)
                while info == Jobs.Pending || info == Jobs.Running || info == Jobs.Submitted
                    sleep(DFC.Service.SLEEP_TIME)
                    info = state(jobdir)
                end
            end
        end
        return true
    end
    run_queue(ARGS[1])
    """)
    
    for (i, f) in enumerate(workflow.steps)
        write(joinpath(qdir, "$i.jl"), @code_string f(job))
    end
end
