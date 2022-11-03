##### JOB INTERACTIONS #########
"""
    load(server::Server, j::Job)

Tries to load the [`Job`](@ref) from `server` at directory `j.dir`.
If no exact matching directory is found, a list of job directories that comprise
`j.dir` will be returned.
"""
function RemoteHPC.load(server::Server, j::Job)
    if j.version == last_version(server, j.dir)
        dir = Jobs.main_job_dir(j.dir)
    elseif !occursin(Jobs.VERSION_DIR_NAME, j.dir) && j.version != -1
        dir = Jobs.version_dir(Jobs.main_job_dir(j.dir), j.version)
    else
        dir = j.dir
    end
    dir = abspath(server, dir)
    
    info = load(server, dir)
    if info isa Vector
        return info
    end
    name = info[2].name
    environment = info[2].environment
    
    remote_calculations = info[2].calculations

    structures = Structure[]
    outcalcs = Calculation[]
    for calc in remote_calculations
        exec = calc.exec
        s = split(calc.args)
        redir_id = findfirst(x -> x == '>', s)
        infile  = redir_id === nothing ? "" : s[redir_id-1]
        outfile = redir_id === nothing ? "" : s[end]
        if Calculations.is_wannier_exec(exec) && !isempty(outcalcs) && outcalcs[end].infile == infile
            Calculations.set_flags!(outcalcs[end], :preprocess => outcalcs[end].run, print=false)
            empty!(outcalcs[end].exec.flags)
        elseif !isempty(infile)
            c = FileIO.calculationparser(exec)(IOBuffer(read(server, joinpath(j, infile))))
            if c.structure !== nothing
                push!(structures, c.structure)
            end
            push!(outcalcs, Calculation(splitext(infile)[1], c.flags, c.data, exec, calc.run, infile, outfile))
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
        pseudos[a] = Structures.Pseudo(server.name, realpath(server, joinpath(dir, "$a")), "") 
    end
    for a in structure.atoms
        a.pseudo = pseudos[a.element.symbol]
    end
    version = j.version == -1 ? last_version(server, j.dir) : j.version
    
    return Job(name, structure, calculations, j.dir, header, version, j.copy_temp_folders, server.name, environment)
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
        resp_job_version = JSON3.read(HTTP.put(server, "/increment_version/" * apath).body, Int)
        job.version = resp_job_version
        @info "Job version: $(curver) => $(resp_job_version)."
    end

    remote_calcs = RemoteHPC.Calculation[]
    for c in job.calculations
        append!(remote_calcs, Calculations.remote_calcs(job, c))
    end
    
    file_buffers = write_calculations(job; kwargs...)
    save(server, job.dir, environment, remote_calcs, name = job.name)
                          
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
        if p.path != linkpath
            if ispath(server, linkpath)
                rm(server, linkpath)
            end
            try
                symlink(server, p.path, linkpath)
            catch
                nothing
            end
        end
    end
    Calculations.rm_tmp_flags!.(job.calculations)
    return job
end

"""
    submit(job::Job)

Saves and launches `job`. 
"""
function RemoteHPC.submit(job::Job, workflow=nothing; kwargs...)
    @assert workflow === nothing "Workflows not implemented yet."
    server = Server(job.server)
    save(job, workflow; kwargs...)
    return submit(server, job.dir)
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
        
#         replacements = fill_execs(server, execs)
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


####### VERSIONING ###################
"""
    versions(job::Job)
    versions(server::Server, jobdir::String)

Returs the valid versions of `job`.
"""
versions(job::Job) = versions(Server(job.server), abspath(job)) 
versions(server::Server, dir::AbstractString) = parse.(Int, readdir(server, joinpath(dir, ".versions")))

"""
    last_version(job::Job)
    last_version(s::Server, jobdir::String)

Returns the last version number of `job`.
"""
last_version(job::Job) = last_version(Server(job.server), abspath(job))
function last_version(server::Server, dir::AbstractString)
    t = versions(server, dir)
    return isempty(t) ? 0 : t[end]
end


"""
    switch_version!(job::Job, version::Int)

Changes the version of `job`.
"""
function switch_version!(job::Job, version::Int)
    allvers = versions(job)
    if !(version in allvers)
        error("Version $version does not exist.")
    end
    tj = load(Server(job.server), Job(Jobs.main_job_dir(job), version=version))
    for f in fieldnames(Job)
        setfield!(job, f, getfield(tj, f))
    end
    return job
end

"""
    rm_version!(job::Job, version::Int)

Removes the `version` of the `job`.
"""
function rm_version!(job::Job, version::Int)
    server = Server(job.server) 
    allvers = versions(job)
    if !(version in allvers)
        error("Version $version does not exist.")
    end

    rm(server, joinpath(job, ".versions", "$version"))    
    if version == job.version
        @warn "Job version is the same as the one to be removed, switching to last known version."
        lv = last_version(job)
        if lv == version
            lv = versions(job)[end-1]
        end
        switch_version!(job, lv)
    end
end

##### EXECS #######
function fill_execs(job::Job, server::Server)
    replacements = fill_execs(server,unique(map(x->x.exec, filter(x->x.run, job.calculations))))
    for (e, rep) in replacements
        for c in job.calculations
            if c.exec == e
                c.exec.dir = rep.dir
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
            n = RemoteHPC.name(server, e)
            if n !== nothing
                e.name = n
            else
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
    for c in calcs
        calculation = job[c]
        p = joinpath(job, calculation.outfile)
        if ispath(server, p)
            if Calculations.isprojwfc(calculation)
                dos_files = filter(x->any(y-> occursin(string(y.name), x), job.structure.atoms), searchdir(job, "pdos_a"))
                extra_files = String[] 
                dir = mkdir(tempname())
                for e in dos_files
                    tfp = joinpath(dir, splitpath(e)[end])
                    write(tfp, read(server, e))
                    push!(extra_files, tfp)
                end
            else
                extra_files = String[]
            end
            try
                out[c] = FileIO.outputdata(calculation, IOBuffer(read(server, p)), extra_files...; extra_parse_funcs = get(extra_parse_funcs, c, Pair{String}[]))
                for strkey in (:initial_structure, :final_structure) 
                    if haskey(out[c], strkey)
                        for a in out[c][strkey].atoms
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

"""
    archive(job::Job, archive_directory::AbstractString, description::String=""; present = nothing, version=job.version)

Archives `job` by copying it's contents to `archive_directory` alongside a `results.jld2` file with all the parseable results as a Dict. `description` will be saved in a `description.txt` file in the `archive_directory`. 
"""
function archive(job::Job, archive_directory::AbstractString, description::String = "";
                 present = nothing)
    @assert !Jobs.isarchived(job) "Job was already archived"
    final_dir = config_path("jobs", "archived", archive_directory)
    @assert !ispath(final_dir) "A archived job already exists in $archive_directory"

    cleanup(job)

    @assert present === nothing "Presenting is currently broken."

    Servers.pull(Server(job.server), abspath(job), final_dir)
    !isempty(description) && write(joinpath(final_dir, "description.txt"), description)
    @info "Archived job at $final_dir. If you're done with this one, it is safe to delete the directory at $(job.dir) on Server(\"$(job.server)\")."
    return nothing
end


# TODO
"""
    cleanup(job::Job)
    
Clean temporary files from the [`Job`](@ref).
"""
function cleanup(job::Job)
    server = Server(job.server)
    td = joinpath(job, Jobs.TEMP_CALC_DIR)
    if ispath(server, td)
        rm(server, joinpath(job, Jobs.TEMP_CALC_DIR))
    end
    for v in versions(job)
        vpath = version_dir(job, v)
        
        if ispath(server, joinpath(vpath, Jobs.TEMP_CALC_DIR))
            rm(server, joinpath(job, Jobs.TEMP_CALC_DIR))
        end
    end
end

