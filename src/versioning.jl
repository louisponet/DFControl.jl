const VERSION_DIR_NAME = ".versions"

function job_versions(dir::String)
    verdir = joinpath(dir, VERSION_DIR_NAME)
    if ispath(verdir)
        return parse.(Int, readdir(joinpath(dir, VERSION_DIR_NAME)))
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

function maybe_increment_version(job::DFJob)
    versions_path = joinpath(job, VERSION_DIR_NAME)
    if !ispath(versions_path)
        mkpath(versions_path)
        job.version = 1
    end
    if ispath(joinpath(job, "job.tt"))
        tjob = DFJob(job.local_dir, version = last_version(job) + 1)
        vpath = version_path(tjob)
        mkpath(vpath)
        cp(job, vpath)
        job.version = last_version(job) + 1
    end
end

"""
    switch_version(job::DFJob, version::Int)

Switches the version of `job` to one of the previously stored ones.
It will save also the current version for future reference.
"""
function switch_version(job::DFJob, version::Int)
    @assert !isrunning(job; print=false) "Can't switch job versions on a running job."
    cur_version = job.version
    if version != cur_version
        verpath = version_path(job.local_dir, version)
        if !ispath(verpath)
            err_str = "Requested job version ($version) is invalid, please choose from:\n"
            for jv in versions(job)
                err_str *= "\t$jv\n"
            end
            @error err_str
        else
            maybe_increment_version(job)
            out = DFJob(verpath)
            curdir = job.local_dir
            mv(out, job.local_dir, force=true)
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

version_assert(job, version) = @assert exists_version(job, version) "Version $version does not exist for job."

"""
    rm_version!(job::DFJob, version::Int)
    rm_versions!(job::DFJob, versions::Int...)

Removes the specified `versions` from the `job` if they exists.
""" 
function rm_version!(job::DFJob, version::Int)
    version_assert(job, version)
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
