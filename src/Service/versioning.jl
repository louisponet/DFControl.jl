using ..DFControl: main_job_dir, VERSION_DIR_NAME
# 0 Job version means that no jobs have ran yet.
"Returns the version found in the .metadata.jld2 if it exists. Otherwise 0."
function main_job_version(dir::AbstractString)
    maindir = main_job_dir(dir)
    mdatapath = joinpath(maindir, ".metadata.jld2")
    if ispath(mdatapath)
        metadata = JLD2.load(mdatapath)
        if haskey(metadata, "version")
            return metadata["version"]
        end
    end
    return 0
end
main_job_version(job::DFJob) = main_job_version(job.local_dir)

function job_versions(dir::AbstractString)
    versions = Int[]
    mainver = main_job_version(dir)
    mainver != 0 && push!(versions, mainver)
    verdir = joinpath(dir, VERSION_DIR_NAME)
    if ispath(verdir)
        append!(versions, parse.(Int, readdir(joinpath(dir, VERSION_DIR_NAME))))
    end
    return sort(unique(versions))
end

"""
    versions(job::DFJob)

Returs the valid versions of `job`.
"""
versions(job::DFJob) = job_versions(main_job_dir(job))
version(job::DFJob) = job.version

function last_job_version(dir::AbstractString)
    versions = job_versions(dir)
    return isempty(versions) ? 0 : versions[end]
end

"""
    last_version(job::DFJob)

Returns the last version number of `job`.
"""
last_version(job::DFJob) = last_job_version(main_job_dir(job))

function version_dir(dir::AbstractString, version::Int)
    tpath = joinpath(dir, VERSION_DIR_NAME, "$version")
    ispath(tpath) && return tpath

    if main_job_version(dir) == version
        return dir
    end
    return error("Version $version not found.")
end
version_dir(job::DFJob) = version_dir(main_job_dir(job), job.version)
version_dir(job::DFJob, version::Int) = version_dir(main_job_dir(job), version)

exists_version(dir::AbstractString, version::Int) = version ∈ job_versions(dir)
exists_version(job::DFJob, version::Int) = exists_version(main_job_dir(job), version)

"""
    maybe_cp_main_version(job::DFJob)

Looks in the `job.local_dir` for the version of the job in the main directory, and copies it to the
respective directory in the `.versions`.
"""
function maybe_cp_main_version(job::DFJob)
    maindir = main_job_dir(job)
    if ispath(joinpath(maindir, "job.tt"))
        tjob = DFJob(maindir)
        cp(tjob, joinpath(tjob, VERSION_DIR_NAME, "$(job.version)"); force = true)
    end
end

"""
    switch_version!(job::DFJob[, version::Int])

Switches the version of `job` to one of the previously stored ones.
It will save also the current version for future reference.
"""
function switch_version!(job::DFJob, version::Int)
    cur_version = job.version
    if version != cur_version
        version_assert(job, version)
        out = DFJob(main_job_dir(job); version = version)
        for f in fieldnames(DFJob)
            setfield!(job, f, getfield(out, f))
        end
    end
    return job
end

function switch_version!(job::DFJob)
    vs = versions(job)
    timestamps = []
    for v in vs
        mdatapath = joinpath(version_dir(job, v), ".metadata.jld2")
        if ispath(mdatapath) &&
           haskey(load(mdatapath), "metadata") &&
           haskey(load(mdatapath)["metadata"], :timestamp)
            push!(timestamps,
                  string(round(load(mdatapath)["metadata"][:timestamp], Dates.Second)))
        else
            push!(timestamps, "")
        end
    end

    menu = RadioMenu(join.(zip(["v$v" for v in string.(vs)], timestamps),
                           ("\tsaved on: ",)))
    choice = request("Please select which version to switch to:", menu)
    if choice != -1
        return switch_version!(job, vs[choice])
    end
end

function version_assert(job, version)
    @assert exists_version(job, version) "Version $version does not exist for job."
end

"""
    rm_version!(job::DFJob, version::Int)
    rm_versions!(job::DFJob, versions::Int...)

Removes the specified `versions` from the `job` if they exist.
"""
function rm_version!(job::DFJob, version::Int)
    version_assert(job, version)
    if version == main_job_version(job)
        for f in readdir(job.local_dir)
            if occursin(".workflow", f) || f == VERSION_DIR_NAME || splitext(f)[2] == ".jl"
                continue
            else
                rm(joinpath(job, f))
            end
        end
        md = main_job_dir(job)
        lv = last_version(job)
        if lv == version
            lv = versions(job)[end-1]
        end
        if lv != 0
            tj = DFJob(md; version = lv)
            cp(tj, md; force = true)
        end

    else
        rm(version_dir(job, version); recursive = true)
    end
    if version == job.version
        @warn "Job version is the same as the one to be removed, switching to last known version."
        lv = last_version(job)
        if lv == version
            lv = versions(job)[end-1]
        end
        switch_version!(job, lv)
    end
end

function rm_versions!(job::DFJob, versions::Int...)
    for v in versions
        rm_version!(job, v)
    end
end

"Removes temporary directories of the specified job versions."
function rm_tmp_dirs!(job, vers = versions(job)...)
    for v in vers
        version_assert(job, v)
        rm(joinpath(version_dir(job, v), TEMP_CALC_DIR); recursive = true)
    end
end
