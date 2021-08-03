function DFControl.Job(dir::String, s="localhost"; version::Int = -1)
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
    job = JSON3.read(resp.body, Job)
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

function save(job::Job)
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
    resp_job = JSON3.read(HTTP.post(server, "/jobs/" * job.dir, [], JSON3.write(job)).body, Job)
    for f in files_to_copy
        push(f, server, joinpath(job.dir, splitdir(f)[end]))
    end
        
    @info "Job version: $(curver) => $(resp_job.version)."
    for f in fieldnames(Job)
        setfield!(job, f, getfield(resp_job, f))
    end
    if haskey(job.metadata, :timestamp)
        job.metadata[:timestamp] = DateTime(job.metadata[:timestamp])
    end
    rm_tmp_flags!(job)
    return job
end

function submit(job::Job)
    server = maybe_start_server(job)
    verify_execs(job, server)
    save(job)
    HTTP.put(server, "/jobs/" * job.dir)
end    

"Runs some checks on the set flags for the calculations in the job, and sets metadata (:prefix, :outdir etc) related flags to the correct ones. It also checks whether flags in the various calculations are allowed and set to the correct types."
function sanitize_flags!(job::Job)
    set_flags!(job, :prefix => "$(job.name)"; print = false)
    if any(x -> eltype(x) == Wannier90, job.calculations) && any(DFC.isnscf, job.calculations)
        nscfcalc = getfirst(DFC.isnscf, job.calculations)
        if eltype(nscfcalc) == Elk
            DFC.set_flags!(job, :num_bands => length(nscfcalc[:wann_bands]); print=false)
            nscfcalc[:wann_projections] = DFC.projections_string.(unique(filter(x -> !isempty(projections(x)), atoms(job))))
            nscfcalc[:elk2wan_tasks]    = ["602", "604"]
            nscfcalc[:wann_seedname]    = Symbol(job.name)
            if job[:wannier_plot] == true
                push!(nscfcalc[:elk2wan_tasks], "605")
            end
        end
    end
    for i in filter(x -> eltype(x) == QE, job.calculations)
        outdir = joinpath(job, DFC.TEMP_CALC_DIR)
        DFC.set_flags!(i, :outdir => "$outdir"; print = false)
    end
    return sanitize_flags!.(job.calculations, (job.structure,))
end


function sanitize_flags!(c::Calculation{QE}, structure::DFC.Structure)
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

    return DFC.convert_flags!(c)
end

function set_hubbard_flags!(c::Calculation{QE}, str::DFC.Structure{T}) where {T}
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

function set_starting_magnetization_flags!(c::Calculation{QE},
                                           str::DFC.Structure{T}) where {T}
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
    sanitize_flags!(calculation::Calculation, str::DFC.Structure)

Cleans up flags, i.e. remove flags that are not allowed and convert all
flags to the correct types.
Tries to correct common errors for different calculation types.
"""
function sanitize_flags!(calculation::Calculation, str::DFC.Structure)
    return convert_flags!(calculation)
end


function rm_tmp_flags!(job::Job)
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

function sanitize_pseudos!(job::Job)
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

function sanitize_magnetization!(job::Job)
    if !any(x -> eltype(x) == QE, job.calculations)
        return
    end
    return DFC.sanitize_magnetization!(job.structure)
end

function find_cutoffs(job::Job)
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

"""
    isrunning(job::Job)

Returns whether a job is running or not. If the job was
submitted using `slurm`, a `QUEUED` status also counts as
running.
"""
function isrunning(job::Job)
    server = maybe_start_server(job)
    return JSON3.read(HTTP.get(server, "/job_isrunning/" * job.dir).body, Bool)
end

"""
    versions(job::Job)

Returs the valid versions of `job`.
"""
function versions(job::Job)
    server = maybe_start_server(job)
    return JSON3.read(HTTP.get(server, "/job_versions/" * job.dir).body, Vector{Int})
end

"""
    last_version(job::Job)

Returns the last version number of `job`.
"""
function last_version(job::Job)
    server = maybe_start_server(job)
    return JSON3.read(HTTP.get(server, "/job_versions/" * job.dir).body, Vector{Int})[end]
end

"""
    last_running_calculation(job::Job)

Returns the last `Calculation` for which an output file was created.
"""
function last_running_calculation(job::Job)
    server = maybe_start_server(job)
    resp = HTTP.get(server, "/last_running_calculation", [], JSON3.write(job))
    if resp.status == 204
        return nothing
    else
        return job[JSON3.read(resp.body, Int)]
    end
end

function outputdata(job::Job)
    server = maybe_start_server(job)
    tmp_path = JSON3.read(HTTP.get(server, "/outputdata", [], JSON3.write(job)).body, String)
    local_temp = tempname() *".jld2"
    pull(server, tmp_path, local_temp)
    return DFC.JLD2.load(local_temp)["outputdata"]
end

function verify_execs(job::Job, server::Server)
    for e in unique(execs(job))
        if !JSON3.read(HTTP.get(server, "/verify_exec/", [], JSON3.write(e)).body, Bool)
            error("$e is not a valid executable on server $(server.name)")
        end
    end
end

"""
    set_flow!(job::Job, should_runs::Pair{String, Bool}...)

Sets whether or not calculations should be scheduled to run.
The `name` of each calculation in the job will be checked against the string in each pair of `should_runs`, and the
`calculation.run` will be set accordingly.

Example:
```julia
set_flow!(job, "" => false, "scf" => true)
```
would un-schedule all calculations in the job, and schedule the "scf" and "nscf" calculations to run.
"""
function set_flow!(job::Job, should_runs...)
    for (name, run) in should_runs
        for calculation in filter(x -> occursin(name, x.name), job.calculations)
            calculation.run = run
        end
    end
    return job
end

"""
    set_headerword!(job::Job, old_new::Pair{String, String})

Replaces the specified word in the header with the new word.
"""
function set_headerword!(job::Job, old_new::Pair{String,String}; print = true)
    for (i, line) in enumerate(job.header)
        if occursin(first(old_new), line)
            job.header[i] = replace(line, old_new)
            s = """Old line:
            $line
            New line:
            $(job.header[i])
            """
            print && (@info s)
        end
    end
    return job
end

"""
    set_dir!(job::Job, dir::AbstractString; copy=false)
    
Sets `job.dir` to `dir`. If necessary the directory will be created upon saving the job.
If `copy` is set to `true`, all previous calculations and output files of the current job version
(i.e. those in the main job directory) will be copied to the new directory, including the
`outputs` directory with temporary files created during jobs runs.
"""
function set_dir!(job::Job, dir::AbstractString; copy = false)
    if !isabspath(dir)
        dir = joinpath(Server(job), dir)
    end
    if dir[end] == '/'
        dir = dir[1:end-1]
    end
    if copy
        error("TODO: Implement for server side copying")
        # mkpath(dir)
        # cp(job, dir; temp = true)
    end
    job.dir = dir
    for i in job.calculations
        set_dir!(i, dir)
    end
    return job
end

#TODO: work this
"""
    abort(job::Job)

Will try to remove the job from the scheduler's queue.
If the last running calculation happened to be a `Calculation{QE}`, the correct abort file will be written.
For other codes the process is not smooth, and restarting is not guaranteed.
"""
function abort(job::Job)
    lastrunning = job.calculations[last_running_calculation(job)]
    if lastrunning == nothing
        error("Is this job running?")
    end
    if package(lastrunning) == QE
        length(filter(x -> eltype(x) == QE, job.calculations)) > 1 &&
            @warn "It's absolutely impossible to guarantee a graceful abort of a multi job script with QE."

        abortpath = writeabortfile(job, lastrunning)
        while ispath(abortpath)
            continue
        end
        qdel(job)
    else
        qdel(job)
    end
end

"Reads throught the pseudo files and tries to figure out the correct cutoffs"
set_cutoffs!(job::Job) = set_cutoffs!.(job.calculations, find_cutoffs(job)...)

"""
    set_wanenergies!(job::Job, nscf::Calculation{QE}, Emin::Real; Epad=5.0)

Will set the energy window limits correctly according to the projections specified in the
structure of the job and the specified Emin. The output of `nscf` will be used to determine the
DOS, and what the size of the frozen window needs to be to fit enough bands inside it,
depending on the projections.
"""
function set_wanenergies!(job::Job, nscf::Calculation, Emin::Real; Epad = 5.0)
    wancalcs = filter(x->eltype(x) == Wannier90, job.calculations)
    @assert length(wancalcs) != 0 "Job ($(job.name)) has no Wannier90 calculations, nothing to do."
    map(x -> set_wanenergies!(x, structure(job), nscf, Emin; Epad = Epad), wancalcs)
    return job
end

"""
    set_wanenergies!(job::Job, nscf::Calculation{QE}, projwfc::Calculation{QE}, threshold::Real; Epad=5.0)

Will set the energy window limits correctly according to the projections specified in the
structure of the job. The output of `projwfc` and the `threshold` will be used to determine
the minimum limit of the frozen energy window such that the interesting DOS of inside it exceeds
the threshold. `nscf` will be used to determine the DOS, and what the upper limit of the frozen window
needs to be to fit enough bands inside it, depending on the projections.
"""
function set_wanenergies!(job::Job, nscf::Calculation, projwfc::Calculation,
                          threshold::Real; Epad = 5.0)
    hasoutput_assert(projwfc)
    @assert isprojwfc(projwfc) "Please specify a valid projwfc calculation."
    @assert isnscf(nscf) "Please specify a valid nscf calculation."
    Emin = Emin_from_projwfc(job.structure, projwfc, threshold)
    return set_wanenergies!(job, nscf, Emin; Epad = Epad)
end

"""
    set_wanenergies!(job::Job, min_window_determinator::Real; kwargs...)

Sets the energy windows of wannier calculations based on the `job`.
When a projwfc calculation is present in the `job`, `min_window_determinator` will be used to
determine the threshold value for including a band in the window based on the projections, otherwise
it will be used as the `Emin` value from which to start counting the number of bands needed for all
projections.
"""
function set_wanenergies!(job::Job, min_window_determinator::Real; kwargs...)
    nscf_calculation = getfirst(isnscf, job.calculations)
    projwfc_calculation = getfirst(isprojwfc, job.calculations)
    if projwfc_calculation === nothing || !hasoutput(projwfc_calculation)
        @info "No projwfc calculation found with valid output, using $min_window_determinator as Emin"
        return set_wanenergies!(job, nscf_calculation, min_window_determinator; kwargs...)
    else
        @info "Valid projwfc output found, using $min_window_determinator as the dos threshold."
        return set_wanenergies!(job, nscf_calculation, projwfc_calculation,
                                min_window_determinator; kwargs...)
    end
end

"""
    bandgap(job::Job, fermi=nothing)

Calculates the bandgap (possibly indirect) around the fermi level.
Uses the first found bands calculation, if there is none it uses the first found nscf calculation.
"""
function bandgap(job::Job, fermi = nothing)
    band_calcs = getfirst.([isbands, isnscf, isscf], (job.calculations,))
    if all(x -> x === nothing, band_calcs)
        error("No valid calculation found to calculate the bandgap.\nMake sure the job has either a valid bands or nscf calculation.")
    end
    if fermi === nothing
        fermi_calcs = getfirst.([isnscf, isscf], (job.calculations,))
        if all(x -> x === nothing, band_calcs)
            error("No valid calculation found to extract the fermi level.\nPlease supply the fermi level manually.")
        end
        fermi = maximum(readfermi.(filter(x -> x !== nothing, fermi_calcs)))
    end

    bands = readbands.(filter(x -> x !== nothing, band_calcs))
    return minimum(bandgap.(bands, fermi))
end

function readfermi(job::Job)
    ins = filter(x -> (isscf(x) || isnscf(x)) && hasoutfile(x), job.calculations)
    @assert isempty(ins) !== nothing "Job does not have a valid scf or nscf output."
    for i in ins
        o = outputdata(i)
        if haskey(o, :fermi)
            return o[:fermi]
        end
    end
    @warn "No output files with fermi level found."
    return 0.0
end

function readbands(job::Job)
    calculation = getfirst(x -> isbands(x) && hasoutfile(x), job.calculations)
    if calculation === nothing
        calculation = getfirst(x -> isnscf(x) && hasoutfile(x), job.calculations)
        if calculation === nothing
            @warn "Job does not have a valid bands output."
            return nothing
        end
        @warn "No bands calculation found, return bands from nscf calculation."
        return readbands(calculation)
    end
    return readbands(calculation)
end


#TODO: only for QE 
"Reads the pdos for a particular atom. Only works for QE."
function pdos(job::Job, atsym::Symbol, filter_word = "")
    projwfc = getfirst(isprojwfc, job.calculations)
    ats = atoms(job, atsym)
    @assert length(ats) > 0 "No atoms found with name $atsym."
    scf = getfirst(isscf, job.calculations)
    magnetic = any(ismagnetic, atoms(job)) || ismagnetic(scf)
    soc = issoc(scf)
    return pdos(projwfc, atsym, magnetic, soc, filter_word)
end

pdos(job::Job, atom::Atom, args...) = pdos(job, atom.name, args...)

function pdos(job::Job, atoms::Vector{Atom} = atoms(job), args...)
    t_energies, t_pdos = pdos(job, atoms[1], args...)
    for i in 2:length(atoms)
        t1, t2 = pdos(job, atoms[i], args...)
        t_pdos .+= t2
    end
    return (energies = t_energies, pdos = t_pdos)
end

"""
    last_submission(job::Job)

If a job was ever submitted, the last submission date is returned.
Otherwise 0 date is returned.
"""
function last_submission(job::Job)
    return get(job.metadata, :timestap, DateTime(0))
end

"""
    set_present!(job::Job, func::Function)
    set_present!(job::Job, func::String)
    set_present!(job::Job, func::Expr)

Sets a function with the call signature `func(job)` which can be later called using the [`@present`](@ref) macro.
"""
function set_present!(job::Job, func::Function)
    try 
        str = loaded_modules_string() * @code_string func(job)
        set_present!(job, str)
    catch
        error("Could not generate the source string for the supplied function.\nIf you are running in a Jupyter notebook, please supply the source code as a string, or expression to set_present!.\nDon't forget to include the right using statements.")
    end
        
end
function set_present!(job::Job, func::AbstractString) 
    open(joinpath(job, ".present.jl"), "w") do f
        write(f, func)
    end
end
function set_present!(job::Job, func::Expr)
    funcstr = string(func)
    if funcstr[1:5] == "begin"
        funcstr = funcstr[findfirst(isequal('\n'), funcstr)+1:findlast(isequal('\n'), funcstr)-1]
    end
    set_present!(job, funcstr)
end

"""
    present(job)

Calls a present function if it was previously saved using [`set_present!`](@ref) or [`archive`](@ref). 
"""
macro present(job)
    return esc(quote
        if ispath(joinpath($job, ".present.jl"))
            t = include(joinpath($job, ".present.jl"))
            DFControl.with_logger(DFControl.MinLevelLogger(DFControl.current_logger(), DFControl.Logging.Error)) do 
                t($job)
            end
        else
            @error "No presentation function defined.\n Please set it with `set_present!`."
        end
    end)
end

"""
    archive(job::Job, archive_directory::AbstractString, description::String=""; present = nothing, version=job.version)

Archives `job` by copying it's contents to `archive_directory` alongside a `results.jld2` file with all the parseable results as a Dict. `description` will be saved in a `description.txt` file in the `archive_directory`. A different job version can be copied using the `version` kwarg, and with the `present` kwarg a function can be specified that can be later called with the [`@present`](@ref) macro.
"""
function archive(job::Job, archive_directory::AbstractString, description::String=""; present = nothing, version=job.version)
    @assert !isarchived(job) "Job was already archived"
    final_dir = config_path(".archived", archive_directory)
    @assert !ispath(final_dir) "A archived job already exists in $archive_directory"
    mkpath(final_dir)

    present !== nothing && set_present!(job, present)
    out = outputdata(job)
    tj = deepcopy(job)
    switch_version!(tj, version)
    cp(tj, final_dir)
    set_dir!(tj, final_dir)

    JLD2.save(joinpath(final_dir, "results.jld2"), "outputdata", out)
    
    !isempty(description) && write(joinpath(final_dir, "description.txt"), description)
    push!(JOB_REGISTRY.archived, tj.dir)
    @info "Archived job at $(tj.dir). If you're done with this one, it is safe to delete the directory at $(job.dir)."
    write_job_registry()
    return nothing
end
