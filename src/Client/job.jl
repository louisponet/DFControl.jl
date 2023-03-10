##### JOB INTERACTIONS #########
"""
    load(server::Server, j::Job)

Tries to load the [`Job`](@ref) from `server` at directory `j.dir`.
If no exact matching directory is found, a list of job directories that comprise
`j.dir` will be returned.
"""
function RemoteHPC.load(server::Server, j::Job)
    j.server = server.name
    
    if j.version == last_version(server, j.dir)
        dir = Jobs.main_job_dir(j.dir)
    elseif !occursin(Jobs.VERSION_DIR_NAME, j.dir) && j.version != 0
        dir = Jobs.version_dir(Jobs.main_job_dir(j.dir), j.version)
    else
        dir = j.dir
    end
    dir = abspath(server, dir)
    
    info = load(server, dir)
    if info isa Vector
        return info
    end
    name = info.name
    environment = info.environment
    
    remote_calculations = info.calculations

    structures = Structure[]
    outcalcs = Calculation[]
    for calc in remote_calculations
        e = calc.exec
        s = split(calc.args)
        redir_id = findfirst(x -> x == ">", s)
        infile  = redir_id === nothing ? "" : s[redir_id-1]
        outfile = redir_id === nothing ? "" : s[end]
        if Calculations.is_wannier_exec(e) && !isempty(outcalcs) && outcalcs[end].infile == infile
            Calculations.set_flags!(outcalcs[end], :preprocess => outcalcs[end].run, print=false)
            empty!(outcalcs[end].exec.flags)
        elseif !isempty(infile)
            c = FileIO.calculationparser(e)(IOBuffer(read(server, joinpath(dir, infile))))
            if c.structure !== nothing
                push!(structures, c.structure)
            end
            push!(outcalcs, Calculation(splitext(infile)[1], c.flags, c.data, e, calc.run, infile, outfile))
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

    calculations = outcalcs

    pseudos = Dict{Symbol, Pseudo}()
    
    for a in unique(map(x->x.element.symbol, structure.atoms))
        if ispath(server, joinpath(dir, "$a" * ".UPF"))
            pseudos[a] = Structures.Pseudo(server.name, realpath(server, joinpath(dir, "$a" * ".UPF")), "")
        else
            pseudos[a] = Structures.Pseudo("", "", "")
        end
    end
    for a in structure.atoms
        a.pseudo = pseudos[a.element.symbol]
    end
    version = last_version(server, j.dir)
    
    return Job(name, structure, calculations, j.dir, version, j.copy_temp_folders, server.name, environment.name)
end

function write_calculations(job::Job; fillexecs=true)
    server = Server(job.server)
    t = Threads.@spawn if fillexecs
        fill_execs(job, server)
    end
     
    Structures.sanitize!(job.structure)
    Calculations.sanitize_flags!(job.calculations, job.structure, job.name,
                                 "./"*Jobs.TEMP_CALC_DIR)
    Jobs.sanitize_cutoffs!(job)
    
    calculations = [c.infile => IOBuffer() for c in job.calculations]
    Threads.@threads for i in eachindex(job.calculations)
        c = job.calculations[i]
        write(calculations[i][2], c, job.structure)
    end

    wcalcs = filter(x -> eltype(x) == Wannier90, job.calculations)
    if !isempty(wcalcs)
        nscf = getfirst(x->Calculations.isnscf(x), job.calculations)
        
        @assert nscf !== nothing "No NSCF found to generate pw2wannier90 from."
        @assert eltype(nscf) == QE "Only QE based Wannier90 jobs are supported."
        
        for c in wcalcs
            pwcalc = FileIO.qe_generate_pw2wancalculation(c, nscf)
            b = IOBuffer()
            write(b, pwcalc, job.structure)
            push!(calculations, pwcalc.infile => b)
        end
    end
    fetch(t)
    return calculations
end

"""
    save(job::Job)

Saves the job's calculations and `job.sh` submission script in `job.dir`.
Some sanity checks will be performed on the validity of flags, execs, pseudopotentials, etc.
The job will also be registered for easy retrieval at a later stage.

If a previous job is present in the job directory (indicated by a valid job script),
it will be copied to the `.versions` sub directory as the previous version of `job`,
and the version of `job` will be incremented. 
"""
function RemoteHPC.save(job::Job, workflow = nothing; versioncheck=true, kwargs...)
    @assert workflow === nothing "Workflows not implemented yet."
    # First we check whether the job is trying to be saved in a archived directory, absolutely not allowed
    @assert !isrunning(job) "Can't save a job in a directory where another is running."

    @assert !Jobs.isarchived(job)
    "Not allowed to save a job in a archived directory, please specify a different directory."

    if isempty(job.name)
        @warn "Job had no name, changed it to: noname"
        job.name = "noname"
    end
    job.name = replace(job.name, " " => "_")

    curver = job.version

    apath = abspath(job)

    server = Server(job.server)
    @assert job.environment != "" "Please set job environment."
    environment = load(server, Environment(job.environment))
    @assert environment isa Environment "Environment with name $(job.environment) not found!"

    if versioncheck
        lastver = last_version(job)
        job.version = lastver
        maybe_cp_main_version(job)
        job.version = lastver + 1
        @info "Job version: $(curver) => $(job.version)."
    end
    job.dir = Jobs.main_job_dir(job)
    remote_calcs = RemoteHPC.Calculation[]
    for c in job.calculations
        append!(remote_calcs, Calculations.remote_calcs(job, c))
    end
    
    file_buffers = write_calculations(job; kwargs...)
    save(server, job.dir, environment, remote_calcs, name = job.name)

    for (n, b) in file_buffers
        write(server, joinpath(job.dir, n), take!(b))
    end
                          
    for c in job.calculations
        if c.run
            ofile = joinpath(job, c.outfile)
            if ispath(server, ofile)
                rm(server, ofile)
            end
        end
    end
    
    # Here we lazily pull pseudos that are not located on the server we're sending stuff to
    pseudos = unique(y->y[1], map(x->(x.element.symbol, x.pseudo), job.structure.atoms))
    for (el, p) in pseudos

        if p == Pseudo()
            error("Pseudo not set for element $el. Use `set_pseudos!`.")
        end
        
        if p.server !== job.server
            p.pseudo = isempty(p.pseudo) && !isempty(p.server) && isalive(Server(p.server)) && ispath(Server(p.server), p.path) ? read(Server(p.server), p.path, String) : p.pseudo
            p.server = job.server
            p.path = joinpath(job, "$el.UPF") 
            for a in job.structure[element(el)]
                a.pseudo = p
            end
        end
        # Pseudo was not present yet
        if !isempty(p.pseudo)
            write(server, p.path, p.pseudo)
            p.pseudo = ""
        end
        
        linkpath = joinpath(job, "$el.UPF")
        if p.path != linkpath || !ispath(server, linkpath)
            if ispath(server, linkpath)
                rm(server, linkpath)
            end
            try
                symlink(server, p.path, linkpath)
            catch
                nothing
            end
        end
        if !ispath(server, p.path)
            error("Pseudo at $(p.path) does not really exist. Use `set_pseudos!` to reset them.")
        end
    end
    Calculations.rm_tmp_flags!.(job.calculations)
    return job
end

"""
    submit(job::Job)

Saves and launches `job`. 
"""
function RemoteHPC.submit(job::Job, workflow=nothing; priority=RemoteHPC.DEFAULT_PRIORITY, kwargs...)
    @assert workflow === nothing "Workflows not implemented yet."
    server = Server(job.server)
    save(job, workflow; kwargs...)
    return submit(server, job.dir, priority)
end

# function submit(jobs::Vector{Job}, run = true)
#     # To verify all execs
#     server_names = Server.(unique(map(x -> x.server, jobs)))
#     buckets = [jobs[findall(x -> x.server == s, jobs)] for s in server_names]
#     outbuckets = [Vector{Job}(undef, length(b)) for b in buckets]
#     Threads.@threads for i in 1:length(server_names)
#         server = server_names[i]
#         bucket = buckets[i]
#         outbucket = outbuckets[i]
#         execs = unique(vcat([map(x->x.exec, j.calculations) for j in bucket]...))
        
#         replacemen RemoteHPC.abort,ts = fill_execs(server, execs)
#         Threads.@threads for ij in 1:length(bucket)
#             job = bucket[ij]
#             for c in job.calculations
#                 for (e, rep) in replacements
#                     if c.exec == e
#                         c.exec.dir = rep.dir
#                         c.exec.modules = rep.modules
#                     end
#                 end
#             end
#             outbucket[ij] = save(job)
#             run && HTTP.put(server, "/jobs/" * abspath(job))
#         end
#     end
#     return outbuckets
# end

"""
    abort(job::Job)

Will try to remove the job from the scheduler's queue.
If the last running calculation happened to be a `Calculation{QE}`, the correct abort file will be written.
For other codes the process is not smooth, and restarting is not guaranteed.
"""
RemoteHPC.abort(job::Job) = abort(Server(job.server), abspath(job))

"""
    state(job::Job)

Returns the state of a job.
"""
RemoteHPC.state(job::Job) = state(Server(job.server), abspath(job))

"""
    isrunning(job::Job)
    isrunning(s::Server, jobdir::String)

Returns whether a job is running or not. If the job was
submitted using `slurm`, a `QUEUED` status also counts as
running.
"""
function isrunning(args...)
    s = state(args...)
    return s == RemoteHPC.Running || s == RemoteHPC.Pending || s == RemoteHPC.Submitted
end

submission_time(job::Job) = mtime(Server(job.server), Jobs.scriptpath(job))

"""
    last_running_calculation(job::Job)

Returns the last `Calculation` for which an output file was created.
"""
function last_running_calculation(job::Job)
    server   = Server(job.server)
    times    = mtime.((server,), map(x->joinpath(job, x.outfile), job.calculations))
    p        = sortperm(times)
    return job[p[end]]
end



##### EXECS #######
function fill_execs(job::Job, server::Server)
    replacements = fill_execs(server,unique(map(x->x.exec, filter(x->x.run, job.calculations))))
    for (e, rep) in replacements
        for c in job.calculations
            if c.exec == e
                c.exec.path = rep.path
                c.exec.modules = rep.modules
                c.exec.parallel = rep.parallel
            end
        end
    end
end

function fill_execs(server::Server, execs::Vector)
    replacements = Dict{Exec, Exec}()
    for e in execs
        if !RemoteHPC.exists(server, e)
            try
                save(server, e)
            catch
                possibilities = load(server, e)
                if length(possibilities) == 1
                    replacements[e] = load(server, Exec(possibilities[1]))
                else
                    error("""
                    Exec(\"$(e.name)\") not found on Server(\"$(server.name)\").
                    Either save it, or choose one of the possible substitutions:
                    $possibilities""")
                end
            end
        else
            replacements[e] = load(server, e)
        end
    end
    return replacements
end

    
"""
    outputdata(job::Job; server = job.server, calcs::Vector{String}=String[])
    
Finds the output files for each of the calculations of a [`Job`](@ref), and groups all the parsed data into a dictionary.
"""
function outputdata(job::Job; calcs=map(x->x.name, job.calculations), extra_parse_funcs=Dict())
    server = Server(job.server)
    out = Dict{String, Dict{Symbol, Any}}()
    calculations = map(x->job[x], calcs)
    
    if !Jobs.runslocal(job)
        tdir = tempname()
        RemoteHPC.pull(job, tdir, infiles=false, calcs=calculations)
    else
        tdir = job.dir
    end
    
    for c in calculations
        of = Calculations.outfiles(c)
        main_file = joinpath(tdir, of[1])
        if ispath(main_file)
            extra_files = String[]
            if length(of) > 1
                for f in of[2:end]
                    append!(extra_files, searchdir(tdir, f))
                end
            end
            try
                out[c.name] = FileIO.outputdata(c, IOBuffer(read(main_file)), extra_files...; extra_parse_funcs = get(extra_parse_funcs, c.name, Pair{String}[]))
                for strkey in (:initial_structure, :final_structure) 
                    if haskey(out[c.name], strkey)
                        for a in out[c.name][strkey].atoms
                            a.pseudo.server = job.server
                        end
                    end
                end
            catch e
                @warn "Error while reading output of calculation: $(c)\n" e
            end
        end
    end
    return out
end

"""
    bandgap(job::Job, fermi=nothing)

Calculates the bandgap (possibly indirect) around the fermi level.
Uses the first found bands calculation, if there is none it uses the first found nscf calculation.
"""
function bandgap(job::Job, fermi = nothing, outdat = outputdata(job))
    band_calcs = filter(!isnothing, getfirst.([Calculations.isbands, Calculations.isnscf, Calculations.isscf], (job.calculations,)))
    if isempty(band_calcs)
        error("No valid calculation found to calculate the bandgap.\nMake sure the job has either a valid bands or nscf calculation.")
    end

    if fermi === nothing
        fermi_calcs = filter(x -> (Calculations.isvcrelax(x) ||
                                   Calculations.isscf(x) ||
                                   Calculations.isnscf(x)), job.calculations)

        if isempty(fermi_calcs)
            error("No valid calculation found to extract the fermi level.\nPlease supply the fermi level manually.")
        end
        fermi = maximum(x->get(get(outdat, x.name, Dict()), :fermi, -Inf), fermi_calcs)
    end

    bands = map(x->outdat[x.name][:bands], filter(x -> haskey(outdat, x.name), band_calcs))
    return minimum(bandgap.(bands, fermi))
end

"""
    readfermi(job::Job, outdat=outputdata(job))

Tries to read the fermi level from a valid [`Calculation`](@ref) inside `job`. 
"""
function readfermi(job::Job, outdat=outputdata(job))
    ins = filter(x -> (Calculations.isvcrelax(x) || Calculations.isscf(x) || Calculations.isnscf(x)), job.calculations)
    @assert isempty(ins) !== nothing "Job does not have a valid scf or nscf output."
    for i in ins
        if haskey(outdat, i.name)
            o = outdat[i.name]
            if haskey(o, :fermi)
                return o[:fermi]
            end
        end
    end
    @warn "No output files with fermi level found."
    return 0.0
end

"""
    readbands(job::Job, outdat=outputdata(job))

Tries to read the bands from a bands calculation that is present in `job`.
"""
function readbands(job::Job, outdat=outputdata(job))
    calc = getfirst(x -> Calculations.isbands(x), job.calculations)
    if calc === nothing || !haskey(outdat, calc.name) || !haskey(outdat[calc.name], :bands)
        calc = getfirst(x -> Calculations.isnscf(x), job.calculations)
        if calc === nothing || !haskey(outdat, calc.name)|| !haskey(outdat[calc.name], :bands)
            @warn "Job does not have a valid bands output."
            return nothing
        end
        @warn "No bands calculation found, return bands from nscf calculation."
        return outdat[calc.name][:bands]
    end
    return outdat[calc.name][:bands]
end

"""
    gencalc_wan(job::Job, min_window_determinator::Real, extra_wan_flags...; kwargs...)

Automates the generation of wannier calculations based on the `job`.
When a projwfc calculation is present in the `job`, `min_window_determinator` will be used to
determine the threshold value for including a band in the window based on the projections, otherwise
it will be used as the `Emin` value from which to start counting the number of bands needed for all
projections.
`extra_wan_flags` can be any extra flags for the Wannier90 calculation such as `write_hr` etc.
"""
function Calculations.gencalc_wan(job::Job, min_window_determinator::Real,
                                  extra_wan_flags...; kwargs...)
    nscf_calc = getfirst(x -> Calculations.isnscf(x), job.calculations)
    nscf_calc === nothing && "Please first run an nscf calculation."
    if get(nscf_calc, :nosym, false) != true
        @info "'nosym' flag was not set in the nscf calculation.
                If this was not intended please set it and rerun the nscf calculation.
                This generally gives errors because of omitted kpoints, needed for pw2wannier90.x"
    end
    projwfc_calc = getfirst(x -> Calculations.isprojwfc(x), job.calculations)
    outdat = outputdata(job)
    if projwfc_calc === nothing || !haskey(outdat, projwfc_calc.name)
        @info "No projwfc calculation found with valid output, using $min_window_determinator as Emin"
        Emin = min_window_determinator
    else
        @info "Valid projwfc output found, using $min_window_determinator as the dos threshold."
        Emin = Calculations.Emin_from_projwfc(job.structure, outdat[projwfc_calc.name][:states], outdat[projwfc_calc.name][:bands], min_window_determinator)
    end
    return Calculations.gencalc_wan(nscf_calc, job.structure, outdat[nscf_calc.name][:bands], Emin, extra_wan_flags...; kwargs...)
end

####### VERSIONING ###################
"""
    versions(job::Job)
    versions(server::Server, jobdir::String)

Returs the valid versions of `job`.
"""
function versions(server::Server, dir::AbstractString)
    vpath = joinpath(dir, ".versions")
    if ispath(server, vpath)
        parse.(Int, readdir(server, vpath))
    else
        return Int[]
    end
end

"""
    last_version(job::Job)
    last_version(s::Server, jobdir::String)

Returns the last version number of `job`.
"""
last_version(job::Job) = last_version(Server(job.server), abspath(job))
function last_version(server::Server, dir::AbstractString)
    t = versions(server, dir)
    return isempty(t) ? 0 : t[end] + 1
end

"""
    versions(job::Job)

Returs the valid versions of `job`.
"""
versions(job::Job) = versions(Server(job.server), Jobs.main_job_dir(job))
version(job::Job) = job.version

function version_dir(dir::AbstractString, version::Int)
    tpath = joinpath(dir, Jobs.VERSION_DIR_NAME, "$version")
    return tpath
end
version_dir(job::Job) = version_dir(Jobs.main_job_dir(job), job.version)
version_dir(job::Job, version::Int) = version_dir(Jobs.main_job_dir(job), version)

"""
    maybe_cp_main_version(job::Job)

Looks in the `job.dir` for the version of the job in the main directory, and copies it to the
respective directory in the `.versions`.
"""
function maybe_cp_main_version(job::Job)
    maindir = Jobs.main_job_dir(job)
    server = Server(job.server)
    if ispath(server, joinpath(maindir, "job.sh"))
        mkpath(server, joinpath(job, Jobs.VERSION_DIR_NAME, "$(job.version)"))
        cp(job, joinpath(job, Jobs.VERSION_DIR_NAME, "$(job.version)"))
    end
end

"""
    switch_version!(job::Job[, version::Int])

Switches the version of `job` to one of the previously stored ones.
It will save also the current version for future reference.
"""
function switch_version!(job::Job, version::Int)
    cur_version = job.version
    if version != cur_version
        version_assert(job, version)
        if version == last_version(job)
            out = load(Server(job.server), Job(Jobs.main_job_dir(job)))
        else
            
            out = load(Server(job.server), Job(joinpath(Jobs.main_job_dir(job), Jobs.VERSION_DIR_NAME, "$(version)")))
        end
        for f in fieldnames(Job)
            if f == :server
                continue
            end
            setfield!(job, f, getfield(out, f))
        end
    end
    job.version = version
    return job
end

function version_assert(job, version)
    @assert version in versions(job) "Version $version does not exist for job."
end

"""
    rm_version!(job::Job, version::Int)
    rm_versions!(job::Job, versions::Int...)

Removes the specified `versions` from the `job` if they exist.
"""
function rm_version!(job::Job, version::Int)
    version_assert(job, version)
    server = Server(job.server)
    if version == last_version(job)
        for f in readdir(job)
            if f == Jobs.VERSION_DIR_NAME
                continue
            else
                rm(server, joinpath(job, f))
            end
        end
        md = Jobs.main_job_dir(job)
        lv = last_version(job)
        if lv == version
            lv = versions(job)[end-1]
        end
        if lv != 0
            real_path = version_dir(md, lv)
            for f in readdir(real_path)
                cp(server, joinpath(real_path, f), joinpath(md, f))
            end
        end
    else
        rm(server, version_dir(job, version))
    end
end
