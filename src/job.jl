function write_job_registry()
    open(config_path("job_registry.txt"), "w") do f
        for j in JOB_REGISTRY
            write(f, "$j\n")
        end
    end
end

function cleanup_job_registry()
    stale_ids = findall(!ispath, JOB_REGISTRY)
    if stale_ids !== nothing
        jobs_to_remove = JOB_REGISTRY[stale_ids] 
        message = "Removing $(length(jobs_to_remove)) stale jobs (job folder removed) from the registry:\n"
        for j in jobs_to_remove
            message *= "\t$j\n"
        end
        @warn message
        deleteat!(JOB_REGISTRY, stale_ids)
    end
    write_job_registry()
end

function maybe_register_job(abspath::String)
    if !ispath(abspath)
        push!(JOB_REGISTRY, abspath)
    else
        jid = findfirst(isequal(abspath), JOB_REGISTRY)
        if jid === nothing
            push!(JOB_REGISTRY, abspath)
        end
    end
    write_job_registry()
end

function last_job_version(dir::String)
    verdir = joinpath(dir, "versions")
    if ispath(verdir)
        versions = readdir(verdir)
        return parse(Int, versions[end])
    else
        return 1
    end
end

version_path(dir::String, version::Int) = joinpath(dir, "versions", "$version")
job_versions(dir::String) = readdir(joinpath(dir, "versions"))

#TODO should we also create a config file for each job with stuff like server etc? and other config things,
#      which if not supplied could contain the default stuff?
"""
Represents a full DFT job with multiple input files and calculations.
"""
@with_kw_noshow mutable struct DFJob
    name         ::String
    structure    ::AbstractStructure
    inputs       ::Vector{DFInput} = DFInput[]
    local_dir    ::String = pwd()
    server       ::String = getdefault_server()
    server_dir   ::String = ""
    header       ::Vector{String} = getdefault_jobheader()
    metadata     ::Dict = Dict()
    version      ::Int = last_job_version(local_dir) 
    function DFJob(name, structure, calculations, local_dir, server, server_dir, header, metadata, version)
        if !isabspath(local_dir)
            local_dir = abspath(local_dir)
        end
        if isempty(structure.name)
            structure.name = split(name, "_")[1]
        end
        out = new(name, structure, calculations, local_dir, server, server_dir, header, metadata, version)
        # TODO add check_version, and if it doesn't match the one that's currently in the main directory,
        # load the previous one from the versions
        !occursin("versions", local_dir) && maybe_register_job(local_dir)
        return out
    end
end

#TODO implement abinit
function DFJob(job_name, structure::AbstractStructure, calculations::Vector{<:DFInput}, common_flags::Pair{Symbol, <:Any}...;
                    pseudoset=:default,
                    pseudospecifier="",
                    job_kwargs...)

    shared_flags = typeof(common_flags) <: Dict ? common_flags : Dict(common_flags...)
    for i in calculations
        i.flags = merge(shared_flags, i.flags)
    end
    out = DFJob(name = job_name, structure = structure, inputs = calculations; job_kwargs...)
    set_pseudos!(out, pseudoset, pseudospecifier, print=false)
    return out
end

DFJob(job_name, ciffile::String, calculations::Vector{<:DFInput}, args...; kwargs...) =
    DFJob(job_name, Structure(ciffile, name = job_name), calculations, args... ; kwargs...)

function DFJob(job::DFJob, flagstoset...; cell_=copy(cell(job)), atoms_=copy(atoms(job)), name=job.name,
                                          server_dir = job.server_dir,
                                          local_dir  = job.local_dir,
                                          pseudoset  = nothing,
                                          pseudospecifier = "")
    newjob = deepcopy(job)

    set_cell!(newjob, cell_)
    if pseudoset === nothing
        pseudoset, specifier = getpseudoset(job.structure.atoms[1])
        specifier = pseudospecifier === nothing ? specifier : pseudospecifier
        set_atoms!(newjob, atoms_, pseudoset = pseudoset, pseudospecifier=specifier)
    else
        set_atoms!(newjob, atoms_, pseudoset = pseudoset, pseudospecifier= pseudospecifier)
    end
    set_serverdir!(newjob, server_dir)
    set_localdir!(newjob, local_dir)
    newjob.name = name

    set_flags!(newjob, flagstoset..., print=false)
    return newjob
end

"""
    DFJob(job_dir::String; job_fuzzy = "job", kwargs...)

Loads and returns a local DFJob.
If `job_dir` is not a valid path the JOB_REGISTRY will be scanned for a job with matching directory.
The kwargs will be passed to the `DFJob` constructor.
"""
function DFJob(job_dir::String; job_fuzzy="job", version = nothing, kwargs...)
    if ispath(job_dir)
        real_path = job_dir
    else
        matching_jobs = filter(x -> occursin(job_dir, x), JOB_REGISTRY)
        if length(matching_jobs) == 1
            real_path = matching_jobs[1]
        else
            menu = RadioMenu(matching_jobs)
            choice = request("Multiple matching jobs were found, choose one:", menu)
            if choice != -1
                real_path = matching_jobs[choice]
            else
                return
            end
        end
    end
    real_version = version === nothing ? last_job_version(real_path) : version
    return DFJob(;merge(merge((local_dir=real_path,version=real_version), read_job_inputs(joinpath(real_path, searchdir(real_path, job_fuzzy)[1]))), kwargs)...)

end        

function switch_version(job::DFJob, version)
    cur_version = job.version
    if version != cur_version
        verpath = version_path(job.local_dir, version)
        if !ispath(verpath)
            err_str = "Requested job version ($version) is invalid, please choose from:\n"
            for jv in job_versions(job.local_dir)
                err_str *= "\t$jv\n"
            end
            @error err_str
        else
            maybe_increment_version(job)
            out = DFJob(verpath)
            curdir = job.local_dir
            mv(out, job.local_dir)
            rm(verpath, recursive=true)
            for f in fieldnames(DFJob)
                setfield!(job, f, getfield(out,f))
            end
            set_localdir!(job, curdir)
            job.version = version
        end
    end
    return job
end

name(job) = job.name
#-------------------BEGINNING GENERAL SECTION-------------#
scriptpath(job::DFJob)       = joinpath(job.local_dir, "job.tt")
starttime(job::DFJob)        = mtime(scriptpath(job))

runslocal(job::DFJob)        = job.server        =="localhost"
structure(job::DFJob)        = job.structure
isQEjob(job::DFJob)          = any(x->package(x) == QE, inputs(job))
iswannierjob(job::DFJob)     = any(x->package(x) == Wannier90, inputs(job)) && any(x->isnscfcalc(x), inputs(job))
getnscfcalc(job::DFJob)      = getfirst(x -> isnscfcalc(x), inputs(job))
cell(job::DFJob)             = cell(structure(job))

input(job::DFJob, n::String) = getfirst(x -> occursin(n, name(x)), inputs(job))
inputs(job::DFJob)           = job.inputs

"""
    inputs(job::DFJob, names::Vector)

Returns an array of the inputs that match the names.
"""
inputs(job::DFJob, names::Vector, fuzzy=true) = fuzzy ? filter(x -> any(occursin.(names, name(x))), inputs(job)) : input.(job, names)
inputs(job::DFJob, n::String, fuzzy=true) = inputs(job, [n], fuzzy)
inputs(job::DFJob, package_::Package) = filter(x->package(x)==package_, inputs(job))
inpath(job::DFJob, n) = inpath(input(job,n))
outpath(job::DFJob, n) = outpath(input(job,n))

"Runs some checks on the set flags for the inputs in the job, and sets metadata (:prefix, :outdir etc) related flags to the correct ones. It also checks whether flags in the various inputs are allowed and set to the correct types."
function sanitizeflags!(job::DFJob)
    set_flags!(job, :prefix => "$(job.name)", print=false)
    if iswannierjob(job)
	    nscfcalc = getnscfcalc(job)
	    if package(nscfcalc) == QE && hasflag(nscfcalc, :nbnd)
	        set_flags!(job, :num_bands => nscfcalc[:nbnd], print=false)
        elseif package(nscfcalc) == Elk
	        setflags!(job, :num_bands => length(nscfcalc[:wann_bands]))
	        nscfcalc[:wann_projections] = projections_string.(unique(filter(x->!isempty(projections(x)), atoms(job))))
	        nscfcalc[:elk2wan_tasks]    = ["602", "604"]
	        nscfcalc[:wann_seedname]    = Symbol(name(job))
	        if job[:wannier_plot] == true
		        push!(nscfcalc[:elk2wan_tasks], "605")
	        end
        end
    end
    for i in filter(x -> package(x) == QE, inputs(job))
        outdir = isempty(job.server_dir) ? joinpath(job, "outputs") : joinpath(job.server_dir, splitdir(job.local_dir)[end], "outputs")
        set_flags!(i, :outdir => "$outdir", print=false)
	    ecutwfc, ecutrho = find_cutoffs(job) # Ideally this should also be at the end stage
	    if exec(i, "pw.x") !== nothing
    	    if !hasflag(i, :ecutwfc)
        	    @info "No energy cutoff was specified in input with name: $(name(i))\nCalculating one from the pseudo files.\nCalculated ecutwfc=$ecutwfc, ecutrho=$ecutrho."
        	    set_flags!(i, :ecutwfc => ecutwfc, :ecuthro => ecutrho, print=false)
    	    end
	    end
    end
    sanitize_pseudos!(job)
    sanitizeflags!.(inputs(job))
end

function sanitize_pseudos!(job::DFJob)
	all_pseudos = pseudo.(atoms(job))
	uni_dirs    = unique(map(x->x.dir, all_pseudos))
	uni_pseudos = unique(all_pseudos)
	pseudo_dir  = length(uni_dirs) == 1 ? uni_dirs[1] : job.local_dir 
	if length(uni_dirs) > 1
		@info "Found pseudos in multiple directories, copying them to job directory"
		for pseudo in uni_pseudos
			cp(path(pseudo), joinpath(job.local_dir, pseudo.name), force=true)
		end
	end
	for p in all_pseudos
		p.dir = pseudo_dir
	end
end

function sanitize_magnetization!(job::DFJob)
    if !any(x->package(x) == QE, inputs(job))
        return
    end
    sanitize_magnetization!(job.structure)
end

function find_cutoffs(job::DFJob)
    @assert job.server == "localhost" "Cutoffs can only be automatically set if the pseudo files live on the local machine."
    pseudofiles = map(x->x.name, filter(!isempty, [pseudo(at) for at in atoms(job)]))
    pseudodirs  = map(x->x.dir, filter(!isempty, [pseudo(at) for at in atoms(job)]))
    @assert !isempty(pseudofiles) "No atoms with pseudo files found."
    @assert !isempty(pseudodirs) "No valid pseudo directories found in the inputs."
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
    if !any(x->!isempty(projections(x)), atoms(job))
        return
    end
    uats = unique(atoms(job))
    projs = unique([name(at) => [p.orb.name for p in projections(at)] for at in uats])
    set_projections!(job, projs...;print=false)
end

"Checks the last created output file for a certain job."
function runninginput(job::DFJob)
    @assert job.server == "localhost" "Intended use for now is locally."
    t = mtime(scriptpath(job))
    for i in reverse(inputs(job))
        p = outpath(i)
        if ispath(p) && mtime(p) > t
            return i
        end
    end
end

"Finds the input corresponding to the name and returns the full output path."
outpath(job::DFJob, n::String) = outpath(input(job,n))

Base.joinpath(job::DFJob, n::AbstractString) = joinpath(job.local_dir, n)

runslocal_assert(job::DFJob) =
    @assert runslocal(job) "This only works if the job runs on `localhost`."

function Base.pop!(job::DFJob, name::String)
    i = findfirst(x -> x.name == name, job.inputs)
    if i === nothing
        error("No input with name $name found.")
    end
    out = job.inputs[i]
    deleteat!(job.inputs, i)
    return out
end

searchdir(job::DFJob, str::AbstractString) = joinpath.((job,), searchdir(job.local_dir, str))

function maybe_increment_version(job::DFJob)
    version_path = joinpath(job, "versions")
    if !ispath(version_path)
        mkpath(version_path)
        return
    else
        if ispath(joinpath(job, "job.tt"))
            vpath = joinpath(version_path, "$(job.version)")
            mkpath(vpath)
            mv(job, vpath)
            job.version = last_job_version(job.local_dir) + 1
        end
    end
end

function Base.mv(job::DFJob, dest::String)
    tjob = DFJob(job.local_dir)
    mv(joinpath(tjob, "job.tt"), joinpath(dest, "job.tt"))
    for i in inputs(tjob)
        mv(i, dest)
    end
    outputs_path = joinpath(job, "outputs")
    if ispath(outputs_path)
        mv(outputs_path, joinpath(dest, "outputs"))
    end
end



