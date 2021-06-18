const TEMP_CALC_DIR = "outputs"
name(job) = job.name
#-------------------BEGINNING GENERAL SECTION-------------#
scriptpath(job::DFJob)       = joinpath(job.local_dir, "job.tt")
starttime(job::DFJob)        = mtime(scriptpath(job))

runslocal(job::DFJob)        = job.server        =="localhost"
structure(job::DFJob)        = job.structure
isQEjob(job::DFJob)          = any(x->package(x) == QE, calculations(job))
iswannierjob(job::DFJob)     = any(x->package(x) == Wannier90, calculations(job)) && any(x->isnscf(x), calculations(job))
getnscfcalc(job::DFJob)      = getfirst(x -> isnscf(x), calculations(job))
cell(job::DFJob)             = cell(structure(job))

calculation(job::DFJob, n::String) = getfirst(x -> occursin(n, name(x)), calculations(job))
calculations(job::DFJob)           = job.calculations

"""
    calculations(job::DFJob, names::Vector)

Returns an array of the calculations that match the names.
"""
calculations(job::DFJob, names::Vector, fuzzy=true) = fuzzy ? filter(x -> any(occursin.(names, name(x))), calculations(job)) : calculation.(job, names)
calculations(job::DFJob, n::String, fuzzy=true) = calculations(job, [n], fuzzy)
calculations(job::DFJob, package_::Package) = filter(x->package(x)==package_, calculations(job))
inpath(job::DFJob, n) = inpath(calculation(job,n))
outpath(job::DFJob, n) = outpath(calculation(job,n))

"Runs some checks on the set flags for the calculations in the job, and sets metadata (:prefix, :outdir etc) related flags to the correct ones. It also checks whether flags in the various calculations are allowed and set to the correct types."
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
    for i in filter(x -> package(x) == QE, calculations(job))
        outdir = isempty(job.server_dir) ? joinpath(job, TEMP_CALC_DIR) : joinpath(job.server_dir, splitdir(job.local_dir)[end], TEMP_CALC_DIR)
        set_flags!(i, :outdir => "$outdir", print=false)
    end
    sanitize_flags!.(calculations(job), (job.structure,))
end

function sanitize_cutoffs!(job)
    # the assumption is that the most important cutoff calculation is the scf/vcrelax that is ran first 
    ψ_cut_calc = getfirst(x -> hasflag(x, ψ_cutoff_flag(x)), calculations(job))
    if ψ_cut_calc !== nothing
        ψcut = ψ_cut_calc[ψ_cutoff_flag(ψ_cut_calc)]
    else
        ψcut, = find_cutoffs(job) # Ideally this should also be at the end stage
        @assert ψcut != 0.0 "No energy cutoff was specified in any calculation, and the calculated cutoff from the pseudopotentials was 0.0.\nPlease manually set one."
        @info "No energy cutoff was specified in the scf calculation.\nCalculated ψcut=$ψcut."
    end
    for i in calculations(job)
        ψflag = ψ_cutoff_flag(i)
        ψflag !== nothing && !hasflag(i, ψflag) && set_flags!(i, ψflag => ψcut, print=false)
    end
    ρ_cut_calc = getfirst(x -> hasflag(x, ρ_cutoff_flag(x)), calculations(job))
    if ρ_cut_calc !== nothing
        ρcut = ρ_cut_calc[ρ_cutoff_flag(ρ_cut_calc)]
        for i in calculations(job)
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
    if !any(x -> package(x) == QE, calculations(job))
        return
    end
    sanitize_magnetization!(job.structure)
end

function find_cutoffs(job::DFJob)
    @assert job.server == "localhost" "Cutoffs can only be automatically set if the pseudo files live on the local machine."
    pseudofiles = map(x->x.name, filter(!isempty, [pseudo(at) for at in atoms(job)]))
    pseudodirs  = map(x->x.dir, filter(!isempty, [pseudo(at) for at in atoms(job)]))
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
    if !any(x->!isempty(projections(x)), atoms(job))
        return
    end
    uats = unique(atoms(job))
    projs = unique([name(at) => [p.orb.name for p in projections(at)] for at in uats])
    set_projections!(job, projs...;print=false)
end

"""
    last_running_calculation(job::DFJob)

Returns the last `DFCalculation` for which an output file was created.
"""
function last_running_calculation(job::DFJob)
    @assert job.server == "localhost" "Intended use for now is locally."
    t = mtime(scriptpath(job))
    for i in reverse(calculations(job))
        p = outpath(i)
        if ispath(p) && mtime(p) > t
            return i
        end
    end
end

"Finds the calculation corresponding to the name and returns the full output path."
outpath(job::DFJob, n::String) = outpath(calculation(job,n))

"""
    joinpath(job::DFJob, args...)

`joinpath(job.local_dir, args...)`.
"""
Base.joinpath(job::DFJob, args...) = joinpath(job.local_dir, args...)

runslocal_assert(job::DFJob) =
    @assert runslocal(job) "This only works if the job runs on `localhost`."

function Base.pop!(job::DFJob, name::String)
    i = findfirst(x -> x.name == name, job.calculations)
    if i === nothing
        error("No calculation with name $name found.")
    end
    out = job.calculations[i]
    deleteat!(job.calculations, i)
    return out
end

searchdir(job::DFJob, str::AbstractString) = joinpath.((job,), searchdir(job.local_dir, str))


for (f, strs) in zip((:cp, :mv), (("copy", "Copies"), ("move", "Moves")))
    @eval begin
    """
        $($f)(job::DFJob, dest::String; all=false, temp=false, kwargs...)

    $($(strs[2])) the contents of `job.local_dir` to `dest`. If `all=true`, it will also $($(strs[1])) the
    `.version` directory with all previous versions. If `temp=true` it will override
    `job.copy_temp_folders` and $($(strs[1])) also the temporary calculation directories.
    The `kwargs...` are passed to `Base.$($f)`.
    """
    function Base.$f(job::DFJob, dest::String; all=false, temp=false, kwargs...)
        for file in readdir(job.local_dir)
            if !all
                if file == VERSION_DIR_NAME
                    continue
                elseif file == TEMP_CALC_DIR && !(temp || job.copy_temp_folders)
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
end

is_slurm_job(job::DFJob) = haskey(job.metadata, :slurmid)

"""
    isrunning(job::DFJob; print=true)

If the job was submitted through a scheduler like `slurm`,
this will return whether the job is queued or running.

**Note:**
For now only `slurm` is supported as scheduler.
"""
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

"""
    filesize(job::DFJob)

Total filesize on disk for a job and all its versions.
"""
Base.filesize(job::DFJob) = dirsize(job.local_dir)

"""
    cleanup(job::DFJob)
    
Cleanup `job.local_dir` interactively.
"""
function cleanup(job::DFJob)
    labels = String[]
    paths = String[]
    for v in versions(job)
        vpath = version_path(job, v)
        s = round(dirsize(vpath)/1e6, digits=3)
        push!(labels, "Version $v:  $s Mb")
        push!(paths, vpath)
        opath = joinpath(vpath, TEMP_CALC_DIR)
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

timestamp(job) = job.metadata[:timestamp]
has_timestamp(job) = haskey(job.metadata, :timestamp)
