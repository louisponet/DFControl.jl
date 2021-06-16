
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
    copy_temp_folders::Bool = false
    function DFJob(name, structure, calculations, local_dir, server, server_dir, header, metadata, version, copy_temp_folders)
        if local_dir[end] == '/'
            local_dir = local_dir[1:end-1]
        end
        if !isabspath(local_dir)
            local_dir = abspath(local_dir)
        end
        if isempty(structure.name)
            structure.name = split(name, "_")[1]
        end
        last_version = last_job_version(local_dir)
        if isempty(metadata)
            mpath = joinpath(local_dir, ".metadata.jld2")
            if ispath(mpath)
                metadata = load(mpath)["metadata"]
            end
        end
        out = new(name, structure, calculations, local_dir, server, server_dir, header, metadata, version, copy_temp_folders)
        # TODO add check_version, and if it doesn't match the one that's currently in the main directory,
        # load the previous one from the versions
        !occursin("versions", local_dir) && maybe_register_job(local_dir)
        return out
    end
end

#TODO implement abinit
function DFJob(job_name, structure::AbstractStructure, calculations::Vector{<:DFInput}, common_flags::Pair{Symbol, <:Any}...;
                    job_kwargs...)

    shared_flags = typeof(common_flags) <: Dict ? common_flags : Dict(common_flags...)
    for i in calculations
        i.flags = merge(shared_flags, i.flags)
    end
    out = DFJob(name = job_name, structure = structure, inputs = calculations; job_kwargs...)
    return out
end

function DFJob(job::DFJob, flagstoset...; cell_=copy(cell(job)), atoms_=copy(atoms(job)), name=job.name,
                                          server_dir = job.server_dir,
                                          local_dir  = job.local_dir)
    newjob = deepcopy(job)

    set_cell!(newjob, cell_)
    set_atoms!(newjob, atoms_)
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
function DFJob(job_dir::String; job_fuzzy="job.tt", version = nothing, kwargs...)
    if !isempty(job_dir) && ispath(abspath(job_dir)) && !isempty(searchdir(abspath(job_dir), job_fuzzy))
        real_path = abspath(job_dir)
    else
        real_path = request_job(job_dir)
        real_path === nothing && return
    end
    real_version = version === nothing ? last_job_version(real_path) : version
    return DFJob(;merge(merge((local_dir=real_path,version=real_version), read_job_inputs(joinpath(real_path, searchdir(real_path, job_fuzzy)[1]))), kwargs)...)

end        

name(job) = job.name
#-------------------BEGINNING GENERAL SECTION-------------#
scriptpath(job::DFJob)       = joinpath(job.local_dir, "job.tt")
starttime(job::DFJob)        = mtime(scriptpath(job))

runslocal(job::DFJob)        = job.server        =="localhost"
structure(job::DFJob)        = job.structure
isQEjob(job::DFJob)          = any(x->package(x) == QE, inputs(job))
iswannierjob(job::DFJob)     = any(x->package(x) == Wannier90, inputs(job)) && any(x->isnscf(x), inputs(job))
getnscfcalc(job::DFJob)      = getfirst(x -> isnscf(x), inputs(job))
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
function sanitize_flags!(job::DFJob)
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
    end
    sanitize_flags!.(inputs(job), (job.structure,))
end

function sanitize_cutoffs!(job)
    # the assumption is that the most important cutoff calculation is the scf/vcrelax that is ran first 
    ψ_cut_calc = getfirst(x -> hasflag(x, ψ_cutoff_flag(x)), inputs(job))
    if ψ_cut_calc !== nothing
        ψcut = ψ_cut_calc[ψ_cutoff_flag(ψ_cut_calc)]
    else
        ψcut, = find_cutoffs(job) # Ideally this should also be at the end stage
        @assert ψcut != 0.0 "No energy cutoff was specified in any input, and the calculated cutoff from the pseudopotentials was 0.0.\nPlease manually set one."
        @info "No energy cutoff was specified in the scf input.\nCalculated ψcut=$ψcut."
    end
    for i in inputs(job)
        ψflag = ψ_cutoff_flag(i)
        ψflag !== nothing && !hasflag(i, ψflag) && set_flags!(i, ψflag => ψcut, print=false)
    end
    ρ_cut_calc = getfirst(x -> hasflag(x, ρ_cutoff_flag(x)), inputs(job))
    if ρ_cut_calc !== nothing
        ρcut = ρ_cut_calc[ρ_cutoff_flag(ρ_cut_calc)]
        for i in inputs(job)
            ρflag = ρ_cutoff_flag(i)
            ρflag !== nothing && set_flags!(i, ρflag => ρcut, print=false)
        end
    end
end

function sanitize_pseudos!(job::DFJob)
    all_pseudos = pseudo.(atoms(job))
    uni_dirs    = unique(map(x->x.dir, all_pseudos))
    uni_pseudos = unique(all_pseudos)
    @assert all(ispath.(path.(uni_pseudos))) "Some Pseudos could not be found, please set them to existing ones."
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
    if !any(x -> package(x) == QE, inputs(job))
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

Base.joinpath(job::DFJob, args...) = joinpath(job.local_dir, args...)

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

for f in (:cp, :mv)
    @eval function Base.$f(job::DFJob, dest::String; all=false, temp=false, kwargs...)
        for file in readdir(job.local_dir)
            if !all
                if file == VERSION_DIR_NAME
                    continue
                elseif file == "outputs" && !(temp || job.copy_temp_folders)
                    continue
                end
            end
            if joinpath(job, file) == abspath(dest)
                continue
            end
            $f(joinpath(job, file), joinpath(dest, file); kwargs...)
        end
    end
end

is_slurm_job(job::DFJob) = haskey(job.metadata, :slurmid)

function isrunning(job::DFJob; print=true)
    if is_slurm_job(job)
        return slurm_isrunning(job)
    end
    print && @warn "Job scheduler unknown."
    return false
end

function dirsize(path::String)
    totsize = 0.0
    for (root, dirs, files) in walkdir(path)
        for file in files
            totsize += filesize(root, file)
        end
    end
    return totsize
end

"Total filesize on disk for a job and all its versions."
Base.filesize(job::DFJob) = dirsize(job.local_dir)

"Cleanup job files interactively."
function cleanup(job::DFJob)
    labels = String[]
    paths = String[]
    for v in versions(job)
        vpath = version_path(job, v)
        s = round(dirsize(vpath)/1e6, digits=3)
        push!(labels, "Version $v:  $s Mb")
        push!(paths, vpath)
        opath = joinpath(vpath, "outputs")
        if ispath(opath)
            s_out = round(dirsize(opath)/1e6, digits=3)
            push!(labels, "Version $v/outputs:  $s_out Mb")
            push!(paths, opath)
        end
    end
    menu = MultiSelectMenu(labels)
    choices = request("Select job files to delete:", menu)
    for i in choices
        if ispath(paths[i]) # Could be that outputs was already deleted
            @info "Deleting $(paths[i])"
            rm(paths[i], recursive=true)
        end
    end
end

save_metadata(job) = jldsave(joinpath(job, ".metadata.jld2"); metadata=job.metadata)
