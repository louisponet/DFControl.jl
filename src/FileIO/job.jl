using ..Servers: scheduler_directive_prefix, scheduler_name_flag

function write_job_header(f, job::Job, environment)
    s = Servers.Server(job.server)
    scheduler_directive = scheduler_directive_prefix(s)
    
    if environment !== nothing
        write(f, "$scheduler_directive --$(scheduler_name_flag(s))=$(job.name) \n")
        for flag in environment.scheduler_flags
            write(f, "$scheduler_directive  $flag\n")
        end
        for flag in environment.exports
            write(f, "export $flag\n")
        end
    end
    job_header = job.header
    for line in job_header
        if occursin("\n", line)
            write(f, line)
        else
            write(f, line * "\n")
        end
    end
    if ispath(s, "/etc/profile")
        write(f, "source /etc/profile\n")
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

infile_outfile_str(c::Calculation) = "< $(c.infile) > $(c.outfile)"

function writetojob(f, job, calculation::Calculation, environment; kwargs...)
    should_run = calculation.run
    if !should_run
        write(f, "#")
    end
    writeexec(f, calculation.exec, environment)
    write(f, infile_outfile_str(calculation))
end

function writetojob(f, job, _calculation::Calculation{Wannier90}, environment; kwargs...)
    filename   = _calculation.infile
    should_run = _calculation.run
    nscf = getfirst(x -> Calculations.isnscf(x), job.calculations)
    
    @assert nscf !== nothing "No NSCF found to generate pw2wannier90 from."
    @assert eltype(nscf) == QE "Only QE based Wannier90 jobs are supported."

    pw2wan_exec = Exec(nscf.exec, exec="pw2wannier90.x")
    empty!(pw2wan_exec.flags)

    preprocess   = get(_calculation, :preprocess, false)

    if !preprocess || !should_run
        write(f, "#")
    end
    writeexec(f, _calculation.exec, environment)
    write(f,
          "-pp $filename > $(_calculation.outfile)\n")

    if !preprocess || !should_run
        write(f, "#")
    end
    writeexec(f, pw2wan_exec, environment)
    write(f, "< pw2wan_$(splitext(filename)[1]).in > pw2wan_$(splitext(filename)[1]).out \n")

    if !should_run
        write(f, "#")
    end
    writeexec(f, _calculation.exec, environment)
    write(f, "$filename > $(_calculation.outfile)")
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
    write(f, job::Job, environment::Environment; kwargs...)

Writes the job script for the job.
Kwargs will be passed down to various writetojob functions.
"""
function Base.write(job_buffer::IO, job::Job, environment::Union{Nothing, Environment}; kwargs...)
    cursize = job_buffer.size

    write(job_buffer, "#!/bin/bash\n")
    write_job_header(job_buffer, job, environment)
    write_job_preamble(job_buffer, job)
    
    # abicalculations = filter(x -> eltype(x) == Abinit, job.calculations)
    # !isempty(abicalculations) && writetojob(job_buffer, job, abicalculations, environment; kwargs...)
    # elkcalculations = filter(x -> eltype(x) == Elk, job.calculations)
    # !isempty(elkcalculations) &&
    #     push!(written_calculations, writetojob(job_buffer, job, elkcalculations, environment; kwargs...)...)

    for i in job.calculations
        writetojob(job_buffer, job, i, environment; kwargs...)
        write(job_buffer, "\n")
    end
    write_job_postamble(job_buffer, job)
    return job_buffer.size - cursize
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
        # if occursin("pw2wannier90", efile)
        #     continue
        # end
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
        qe_parse_calculation
    elseif Calculations.is_wannier_exec(exec)
        wan_parse_calculation
    elseif Calculations.is_elk_exec(exec)
        elk_parse_calculation
    elseif Calculations.is_julia_exec(exec)
        julia_parse_calculation
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
            #TODO : Bad 
            if occursin("#SBATCH", line) || occursin("#HQ", line)
                if occursin("-J", line) || occursin("name=", line)
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
    ues = unique(map(x->x.exec, calcs))
    for e in ues
        t = Database.replacements(e)
        if !isempty(t)
            b = getfirst(x -> abspath(e.dir) == abspath(x.dir) && e.exec == x.exec, t)
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
    deleteat!(header, findall(x -> occursin("source /etc/profile", x), header))
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
        elseif !isempty(calc[:contents])
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

    to_rm = Int[]
    for (ic, c) in enumerate(outcalcs)
        id2 = findnext(x->x.name == c.name, outcalcs, ic+2)
        if id2 !== nothing && Calculations.is_wannier_exec(c.exec)
            push!(to_rm, ic + 1)
            push!(to_rm, id2)
            Calculations.set_flags!(c, :preprocess => c.run, print=false)
            Calculations.set_flags!(c, :wannier_plot => get(outcalcs[ic+1], :write_unk, false), print=false)
        end
    end
    deleteat!(outcalcs, to_rm)

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
