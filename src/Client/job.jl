function DFControl.DFJob(dir::String, s="localhost"; version::Int = -1)
    server = maybe_start_server(s)
    # server = Server(s)
    # dir = dir[1] == '/' ? dir[2:end] : dir
    resp = HTTP.get(server, "/jobs/" * dir, [], JSON3.write(version))
    # Supplied dir was not a valid path, so we ask
    # previously registered jobs on the server that
    # contain dir.
    if resp.status == 204
        dir = request_job_dir(dir, server)
        dir === nothing && return
        resp = HTTP.get(server, "/jobs/" * dir, [], JSON3.write(version))
    end 
    # @show String(resp.body)
    job = JSON3.read(resp.body, DFJob)
    job.server = server.name
    if haskey(job.metadata, :timestamp)
        job.metadata[:timestamp] = DateTime(job.metadata[:timestamp])
    end
    return job
end

function request_job_dir(dir::String, server::Server)
    resp = HTTP.get(server, "/registered_jobs/" *  dir)
    matching_jobs = reverse(JSON3.read(resp.body, Vector{Tuple{String, DateTime}}))
    if length(matching_jobs) == 1
        return matching_jobs[1][1]
    elseif length(matching_jobs) == 0
        error("No jobs found matching $dir")
    elseif isdefined(Base, :active_repl)
        choices = ["$j -- $t" for (j, t) in matching_jobs]
        menu = RadioMenu(choices)
        choice = request("Multiple matching jobs were found, choose one:", menu)
        if choice != -1
            return matching_jobs[choice][1]
        else
            return nothing
        end
    else
        err_msg = "No concrete job for $dir found, closest matches are:"
        for m in choices
            err_msg = join([err_msg, m], "\n")
        end
        error(err_msg)
    end
end

function save(job::DFJob)
    # First we check whether the job is trying to be saved in a archived directory, absolutely not allowed
    @assert !DFC.isarchived(job)
        "Not allowed to save a job in a archived directory, please specify a different directory with `set_dir!"
    server = maybe_start_server(job)
    @assert !isrunning(job) "Can't save a job in a directory where another is running."

    if isempty(job.name)
        @warn "Job had no name, changed it to: noname"
        job.name = "noname"
    end
    
    sanitize_cutoffs!(job)
    sanitize_pseudos!(job)
    sanitize_magnetization!(job)
    sanitize_projections!(job)
    sanitize_flags!(job)
    
    curver = job.version
    resp_job = JSON3.read(HTTP.post(server, "/jobs/" * job.dir, [], JSON3.write(job)).body, DFJob)
    @info "Job version: $(curver) => $(resp_job.version)."
    for f in fieldnames(DFJob)
        setfield!(job, f, getfield(resp_job, f))
    end
    if haskey(job.metadata, :timestamp)
        job.metadata[:timestamp] = DateTime(job.metadata[:timestamp])
    end
    rm_tmp_flags!(job)
    return job
end

function submit(job::DFJob)
    server = maybe_start_server(job)
    DFC.verify_execs(job, server)
    save(job)
end    

"Runs some checks on the set flags for the calculations in the job, and sets metadata (:prefix, :outdir etc) related flags to the correct ones. It also checks whether flags in the various calculations are allowed and set to the correct types."
function sanitize_flags!(job::DFJob)
    set_flags!(job, :prefix => "$(job.name)"; print = false)
    if DFC.iswannierjob(job)
        nscfcalc = DFC.getnscfcalc(job)
        if DFC.package(nscfcalc) == Elk
            set_flags!(job, :num_bands => length(nscfcalc[:wann_bands]))
            nscfcalc[:wann_projections] = DFC.projections_string.(unique(filter(x -> !isempty(projections(x)), atoms(job))))
            nscfcalc[:elk2wan_tasks]    = ["602", "604"]
            nscfcalc[:wann_seedname]    = Symbol(name(job))
            if job[:wannier_plot] == true
                push!(nscfcalc[:elk2wan_tasks], "605")
            end
        end
    end
    for i in filter(x -> DFC.package(x) == QE, job.calculations)
        outdir = joinpath(job, DFC.TEMP_CALC_DIR)
        set_flags!(i, :outdir => "$outdir"; print = false)
    end
    return sanitize_flags!.(job.calculations, (job.structure,))
end

function rm_tmp_flags!(job::DFJob)
    rm_flags!(job, :prefix, :outdir; print=false)
    rm_flags!(job, :nspin; print=false)
end

function sanitize_cutoffs!(job)
    # the assumption is that the most important cutoff calculation is the scf/vcrelax that is ran first 
    ψ_cut_calc = getfirst(x -> hasflag(x, DFC.ψ_cutoff_flag(x)), job.calculations)
    if ψ_cut_calc !== nothing
        ψcut = ψ_cut_calc[DFC.ψ_cutoff_flag(ψ_cut_calc)]
    else
        ψcut, = find_cutoffs(job) # Ideally this should also be at the end stage
        @assert ψcut != 0.0 "No energy cutoff was specified in any calculation, and the calculated cutoff from the pseudopotentials was 0.0.\nPlease manually set one."
        @info "No energy cutoff was specified in the scf calculation.\nCalculated ψcut=$ψcut."
    end
    for i in job.calculations
        ψflag = DFC.ψ_cutoff_flag(i)
        ψflag !== nothing &&
            !hasflag(i, ψflag) &&
            set_flags!(i, ψflag => ψcut; print = false)
    end
    ρ_cut_calc = getfirst(x -> hasflag(x, DFC.ρ_cutoff_flag(x)), job.calculations)
    if ρ_cut_calc !== nothing
        ρcut = ρ_cut_calc[DFC.ρ_cutoff_flag(ρ_cut_calc)]
        for i in job.calculations
            ρflag = DFC.ρ_cutoff_flag(i)
            ρflag !== nothing && set_flags!(i, ρflag => ρcut; print = false)
        end
    end
end

function sanitize_pseudos!(job::DFJob)
    all_pseudos = DFC.pseudo.(atoms(job))
    uni_dirs    = unique(map(x -> x.dir, all_pseudos))
    uni_pseudos = unique(all_pseudos)
    if !all(ispath.(DFC.path.(uni_pseudos)))
        @warn "Some Pseudos could not be located on disk."
    end
    pseudo_dir = length(uni_dirs) == 1 ? uni_dirs[1] : job.dir
    if length(uni_dirs) > 1
        @info "Found pseudos in multiple directories, copying them to job directory"
        for pseudo in uni_pseudos
            cp(DFC.path(pseudo), joinpath(job.dir, pseudo.name); force = true)
        end
    end
    for p in all_pseudos
        p.dir = pseudo_dir
    end
end

function sanitize_magnetization!(job::DFJob)
    if !any(x -> package(x) == QE, job.calculations)
        return
    end
    return DFC.sanitize_magnetization!(job.structure)
end

function find_cutoffs(job::DFJob)
    @assert job.server == "localhost" "Cutoffs can only be automatically set if the pseudo files live on the local machine."
    pseudofiles = map(x -> x.name, filter(!isempty, [DFC.pseudo(at) for at in atoms(job)]))
    pseudodirs  = map(x -> x.dir, filter(!isempty, [DFC.pseudo(at) for at in atoms(job)]))
    @assert !isempty(pseudofiles) "No atoms with pseudo files found."
    @assert !isempty(pseudodirs) "No valid pseudo directories found in the calculations."
    maxecutwfc = 0.0
    maxecutrho = 0.0
    for d in pseudodirs
        for f in pseudofiles
            pth = joinpath(d, f)
            if ispath(pth)
                ecutwfc, ecutrho = read_cutoffs_from_pseudofile(pth)
                if ecutwfc != nothing && ecutrho != nothing
                    maxecutwfc = ecutwfc > maxecutwfc ? ecutwfc : maxecutwfc
                    maxecutrho = ecutrho > maxecutrho ? ecutrho : maxecutrho
                end
            end
        end
    end
    return maxecutwfc, maxecutrho
end

function sanitize_projections!(job::DFJob)
    if !any(x -> !isempty(projections(x)), atoms(job))
        return
    end
    uats = unique(atoms(job))
    projs = unique([name(at) => [p.orb.name for p in projections(at)] for at in uats])
    return set_projections!(job, projs...; print = false)
end


"""
    isrunning(job::DFJob)

Returns whether a job is running or not. If the job was
submitted using `slurm`, a `QUEUED` status also counts as
running.
"""
function isrunning(job::DFJob)
    server = maybe_start_server(job)
    return JSON3.read(HTTP.get(server, "/job_isrunning/" * job.dir).body, Bool)
end

"""
    versions(job::DFJob)

Returs the valid versions of `job`.
"""
function versions(job::DFJob)
    server = maybe_start_server(job)
    return JSON3.read(HTTP.get(server, "/job_versions/" * job.dir).body, Vector{Int})
end

"""
    last_version(job::DFJob)

Returns the last version number of `job`.
"""
function last_version(job::DFJob)
    server = maybe_start_server(job)
    return JSON3.read(HTTP.get(server, "/job_versions/" * job.dir).body, Vector{Int})[end]
end

"""
    last_running_calculation(job::DFJob)

Returns the last `DFCalculation` for which an output file was created.
"""
function last_running_calculation(job::DFJob)
    server = maybe_start_server(job)
    resp = HTTP.get(server, "/last_running_calculation", [], JSON3.write(job))
    if resp.status == 204
        return nothing
    else
        return job[JSON3.read(resp.body, Int)]
    end
end

function outputdata(job::DFJob)
    server = maybe_start_server(job)
    tmp_path = JSON3.read(HTTP.get(server, "/outputdata", [], JSON3.write(job)).body, String)
    local_temp = tempname() *".jld2"
    pull(server, tmp_path, local_temp)
    return DFC.JLD2.load(local_temp)["outputdata"]
end

"""
    pull(server::Server, server_file::String, local_file::String)

Pulls `server_file` from the server the `local_file`.
"""
function pull(server::Server, server_file::String, filename::String)
    if server.name == "localhost"
        cp(server_file, filename, force=true)
    else
        run(`scp $(ssh_string(server) * ":" * server_file) $filename`)
    end
end
