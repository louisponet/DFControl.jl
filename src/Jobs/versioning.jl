# 0 Job version means that no jobs have ran yet.
"Returns the version found in the .metadata.jld2 if it exists. Otherwise 0."
main_job_version(dir::AbstractString) = job_versions(dir)[end]
main_job_version(job::Job) = main_job_version(abspath(job))

function job_versions(dir::AbstractString)
    verdir = joinpath(dir, VERSION_DIR_NAME)
    if ispath(verdir)
        versions = parse.(Int, readdir(verdir))
        push!(versions, versions[end]+1)
    else
        versions = Int[0]
    end
       
    return sort(unique(versions))
end

"""
    versions(job::Job)

Returs the valid versions of `job`.
"""
versions(job::Job) = job_versions(main_job_dir(job))
version(job::Job) = job.version
version(dir::AbstractString) = occursin(VERSION_DIR_NAME, dir) ? parse(Int, splitpath(dir)[end]) : main_job_version(dir)

function last_job_version(dir::AbstractString)
    versions = job_versions(dir)
    return isempty(versions) ? 0 : versions[end]
end

"""
    last_version(job::Job)

Returns the last version number of `job`.
"""
last_version(job::Job) = last_job_version(main_job_dir(job))

function version_dir(dir::AbstractString, version::Int)
    tpath = joinpath(dir, VERSION_DIR_NAME, "$version")
    return tpath
end
version_dir(job::Job) = version_dir(main_job_dir(job), job.version)
version_dir(job::Job, version::Int) = version_dir(main_job_dir(job), version)

exists_version(dir::AbstractString, version::Int) = version ∈ job_versions(dir)
exists_version(job::Job, version::Int) = exists_version(main_job_dir(job), version)

"""
    maybe_cp_main_version(job::Job)

Looks in the `job.dir` for the version of the job in the main directory, and copies it to the
respective directory in the `.versions`.
"""
function maybe_cp_main_version(job::Job)
    maindir = main_job_dir(job)
    if ispath(joinpath(maindir, "job.sh"))
        tjob = load_job(maindir)
        cp(tjob, joinpath(tjob, VERSION_DIR_NAME, "$(job.version)"); force = true)
    end
end

"""
    switch_version!(job::Job[, version::Int])

Switches the version of `job` to one of the previously stored ones.
It will save also the current version for future reference.
"""
function switch_version!(job::Job, version::Int)
    cur_version = job.version
    if version != cur_version
        version_assert(job, version)
        out = load_job(main_job_dir(job); version = version)
        for f in fieldnames(Job)
            if f == :server
                continue
            end
            setfield!(job, f, getfield(out, f))
        end
    end
    return job
end

function switch_version!(job::Job)
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
    rm_version!(job::Job, version::Int)
    rm_versions!(job::Job, versions::Int...)

Removes the specified `versions` from the `job` if they exist.
"""
function rm_version!(job::Job, version::Int)
    version_assert(job, version)
    if version == main_job_version(job)
        for f in readdir(job)
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
            real_path = version_dir(md, lv)
            for f in readdir(real_path)
                cp(joinpath(real_path, f), joinpath(md, f), force=true)
            end
        end
    else
        rm(version_dir(job, version); recursive = true)
    end
end

"Removes temporary directories of the specified job versions."
function rm_tmp_dirs!(job, vers = versions(job)...)
    for v in vers
        version_assert(job, v)
        rm(joinpath(version_dir(job, v), TEMP_CALC_DIR); recursive = true)
    end
end
