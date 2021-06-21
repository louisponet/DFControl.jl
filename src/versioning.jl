const VERSION_DIR_NAME = ".versions"

function job_versions(dir::String)
    verdir = joinpath(dir, VERSION_DIR_NAME)
    if ispath(verdir)
        return sort(parse.(Int, readdir(joinpath(dir, VERSION_DIR_NAME))))
    else
        return Int[]
    end
end

"""
    versions(job::DFJob)

Returs the valid versions of `job`.
"""
versions(job::DFJob) = job_versions(job.local_dir)
version(job::DFJob) = job.version

function last_job_version(dir::String)
    versions = job_versions(dir)
    return isempty(versions) ? 0 : versions[end]
end
last_version(job::DFJob) = last_job_version(job.local_dir)

version_path(dir::String, version::Int) = joinpath(dir, VERSION_DIR_NAME, "$version")
version_path(job::DFJob) = version_path(job.local_dir, job.version)
version_path(job::DFJob, version::Int) = version_path(job.local_dir, version)

exists_version(dir::String, version::Int) = version ∈ job_versions(dir)
exists_version(job::DFJob, version::Int) = exists_version(job.local_dir, version)

"""
    maybe_cp_prev_version(job::DFJob)

Looks in the `job.local_dir` for a previous version of the job, and copies it to the
respective directory in the `.versions`.
"""
function maybe_cp_prev_version(job::DFJob)
    versions_path = joinpath(job, VERSION_DIR_NAME)
    if !ispath(versions_path)
        mkpath(versions_path)
    end
    if ispath(joinpath(job, "job.tt"))
        tjob = DFJob(job.local_dir)
        vpath = version_path(tjob)
        cp(tjob, vpath, force=true)
    end
end

"""
    switch_version(job::DFJob[, version::Int])

Switches the version of `job` to one of the previously stored ones.
It will save also the current version for future reference.
"""
function switch_version(job::DFJob, version::Int)
    @assert !isrunning(job; print=false) "Can't switch job versions on a running job."
    cur_version = job.version
    if version != cur_version
        verpath = version_path(job.local_dir, version)
        version_assert(job, version)
        maybe_cp_prev_version(job)
        clean_local_dir!(job)
        out = DFJob(verpath)
        curdir = job.local_dir
        cp(out, job.local_dir, force=true)
        for f in fieldnames(DFJob)
            setfield!(job, f, getfield(out,f))
        end
        set_localdir!(job, curdir)
        job.version = version
    end
    return job
end

function switch_version(job::DFJob)
    vs = versions(job)
    timestamps = []
    for v in vs
        mdatapath = joinpath(version_path(job, v), ".metadata.jld2") 
        if ispath(mdatapath) && haskey(load(mdatapath), "metadata") && haskey(load(mdatapath)["metadata"], :timestamp)
            push!(timestamps, string(load(mdatapath)["metadata"][:timestamp]))
        else
            push!(timestamps, "")
        end
    end
            
        
    menu = RadioMenu(join.(zip(string.(vs), timestamps), (" ",)))
    choice = request("Please select which version to switch to:", menu)
    if choice != -1
        return switch_version(job, vs[choice])
    end
end

version_assert(job, version) = @assert exists_version(job, version) "Version $version does not exist for job."

"""
    rm_version!(job::DFJob, version::Int)
    rm_versions!(job::DFJob, versions::Int...)

Removes the specified `versions` from the `job` if they exist.
""" 
function rm_version!(job::DFJob, version::Int)
    version_assert(job, version)
    if version == job.version
        @warn "Job version is the same as the one to be removed, switching to last known version."
        lv = last_version(job)
        if lv == version
            lv = versions(job)[end-1]
        end
        switch_version(job, lv)
    end
    rm(version_path(job, version), recursive=true)
end

function rm_versions!(job::DFJob, versions::Int...)
    for v in versions
        rm_version!(job, v)
    end
end

"Removes temporary directories of the specified job versions."
function rm_tmp_dirs!(job, vers=versions(job)...)
    for v in vers
        version_assert(job, v)
        rm(joinpath(version_path(job, v), TEMP_CALC_DIR), recursive=true)
    end
end
