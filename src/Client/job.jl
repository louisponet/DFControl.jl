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
    files_to_copy = sanitize_pseudos!(job)
    sanitize_magnetization!(job)
    sanitize_projections!(job)
    sanitize_flags!(job)
    
    curver = job.version
    resp_job = JSON3.read(HTTP.post(server, "/jobs/" * job.dir, [], JSON3.write(job)).body, DFJob)
    for f in files_to_copy
        push(f, server, joinpath(job.dir, splitdir(f)[end]))
    end
        
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
    verify_execs(job, server)
    save(job)
    HTTP.put(server, "/jobs/" * job.dir)
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

"Runs through all the set flags and checks if they are allowed and set to the correct value"
function convert_flags!(calculation::DFCalculation)
    for (flag, value) in DFC.flags(calculation)
        flagtype_ = DFC.flagtype(calculation, flag)
        if flagtype_ == Nothing
            @warn "Flag $flag was not found in allowed flags for exec $(execs(calculation)[2]). Removing flag."
            rm_flags!(calculation, flag)
            continue
        end
        if !(isa(value, flagtype_) || eltype(value) <: flagtype_)
            try
                if isbitstype(eltype(value))
                    if length(value) > 1
                        calculation[flag] = convert(flagtype_, value)
                    else
                        calculation[flag] = convert(eltype(flagtype_), value)
                    end
                else
                    calculation[flag] = convert.(flagtype_, value)
                end
            catch
                error("Input $(name(calculation)): Could not convert :$flag of value $value to the correct type ($flagtype_), please set it to the correct type.")
            end
        end
    end
end

function sanitize_flags!(c::DFCalculation{QE}, structure::DFC.AbstractStructure)
    if DFC.isvcrelax(c)
        #this is to make sure &ions and &cell are there in the calculation 
        !hasflag(c, :ion_dynamics) && set_flags!(c, :ion_dynamics => "bfgs"; print = false)
        !hasflag(c, :cell_dynamics) &&
            set_flags!(c, :cell_dynamics => "bfgs"; print = false)
    end
    #TODO add all the required flags
    if exec(c, "pw.x") !== nothing
        @assert hasflag(c, :calculation) "Please set the flag for calculation with name: $(name(c))"
    end
    # setting hubbard and magnetization flags
    set_hubbard_flags!(c, structure)
    set_starting_magnetization_flags!(c, structure)

    # setting hubbard flags 
    pseudo_dir = DFC.pseudo(atoms(structure)[1]).dir # Pseudos should be all sanitized by now
    set_flags!(c, :pseudo_dir => pseudo_dir; print = false)

    return convert_flags!(c)
end

function set_hubbard_flags!(c::DFCalculation{QE}, str::DFC.AbstractStructure{T}) where {T}
    u_ats = unique(atoms(str))
    isdftucalc = any(x -> dftu(x).U != 0 ||
                              dftu(x).J0 != 0.0 ||
                              sum(dftu(x).J) != 0 ||
                              sum(dftu(x).α) != 0, u_ats) || hasflag(c, :Hubbard_parameters)
    isnc = DFC.isnoncolin(str)
    if isdftucalc
        Jmap = map(x -> copy(dftu(x).J), u_ats)
        Jdim = maximum(length.(Jmap))
        Jarr = zeros(Jdim, length(u_ats))
        for (i, J) in enumerate(Jmap)
            diff = Jdim - length(J)
            if diff > 0
                for d in 1:diff
                    push!(J, zero(eltype(J)))
                end
            end
            Jarr[:, i] .= J
        end
        set_flags!(c, :lda_plus_u    => true, :Hubbard_U     => map(x -> dftu(x).U, u_ats),
                   :Hubbard_alpha => map(x -> dftu(x).α, u_ats),
                   :Hubbard_beta  => map(x -> dftu(x).β, u_ats), :Hubbard_J     => Jarr,
                   :Hubbard_J0    => map(x -> dftu(x).J0, u_ats); print = false)
        isnc && set_flags!(c, :lda_plus_u_kind => 1; print = false)
    else
        rm_flags!(c, :lda_plus_u, :lda_plus_u_kind, :Hubbard_U, :Hubbard_alpha,
                  :Hubbard_beta, :Hubbard_J, :Hubbard_J0, :U_projection_type; print = false)
    end
end

function set_starting_magnetization_flags!(c::DFCalculation{QE},
                                           str::DFC.AbstractStructure{T}) where {T}
    u_ats = unique(atoms(str))
    mags = magnetization.(u_ats)
    starts = T[]
    θs = T[]
    ϕs = T[]
    ismagcalc = DFC.ismagnetic(str)
    isnc = DFC.isnoncolin(str)
    if (ismagcalc && isnc) || (flag(c, :noncolin) !== nothing && flag(c, :noncolin))
        for m in mags
            tm = normalize(m)
            if norm(m) == 0
                push!.((starts, θs, ϕs), 0.0)
            else
                θ = acos(tm[3]) * 180 / π
                ϕ = atan(tm[2], tm[1]) * 180 / π
                start = norm(m)
                push!(θs, θ)
                push!(ϕs, ϕ)
                push!(starts, start)
            end
        end
        set_flags!(c, :noncolin => true; print = false)
        rm_flags!(c, :nspin; print = false)
    elseif ismagcalc
        for m in mags
            push!.((θs, ϕs), 0.0)
            if norm(m) == 0
                push!(starts, 0)
            else
                push!(starts, sign(sum(m)) * norm(m))
            end
        end
        set_flags!(c, :nspin => 2; print = false)
    end
    return set_flags!(c, :starting_magnetization => starts, :angle1 => θs, :angle2 => ϕs;
                      print = false)
end


#TODO implement abinit and wannier90
"""
    sanitize_flags!(calculation::DFCalculation, str::DFC.AbstractStructure)

Cleans up flags, i.e. remove flags that are not allowed and convert all
flags to the correct types.
Tries to correct common errors for different calculation types.
"""
function sanitize_flags!(calculation::DFCalculation, str::DFC.AbstractStructure)
    return convert_flags!(calculation)
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
    s = Server(job)
    pseudo_paths = DFC.path.(uni_pseudos)
    if !all(x -> JSON3.read(HTTP.get(s, "/get_ispath/" * x).body, Bool), pseudo_paths)
        if all(ispath, pseudo_paths)
            @info "Some Pseudos could not be found on Server $(s.name), pushing them from local storage to job dir."
            for p in all_pseudos
                p.dir = job.dir
            end
            return pseudo_paths
        else
            @warn "Some pseudos could not be found locally, and neither on the remote server."
        end
    end
    return String[]
end

function sanitize_magnetization!(job::DFJob)
    if !any(x -> DFC.package(x) == QE, job.calculations)
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

function verify_execs(job::DFJob, server::Server)
    for e in unique(execs(job))
        if !JSON3.read(HTTP.get(server, "/verify_exec/", [], JSON3.write(e)).body, Bool)
            error("$e is not a valid executable on server $(server.name)")
        end
    end
end

