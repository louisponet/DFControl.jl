const VERSION_DIR_NAME = ".versions"
const TEMP_CALC_DIR = "outputs"

name(job) = job.name
#-------------------BEGINNING GENERAL SECTION-------------#
scriptpath(job::DFJob) = joinpath(job.dir, "job.tt")
starttime(job::DFJob)  = mtime(scriptpath(job))

runslocal(job::DFJob)    = job.server == "localhost"
structure(job::DFJob)    = job.structure
isQEjob(job::DFJob)      = any(x -> package(x) == QE, job.calculations)
iswannierjob(job::DFJob) = any(x -> package(x) == Wannier90, job.calculations) && any(x -> isnscf(x), job.calculations)
getnscfcalc(job::DFJob)  = getfirst(x -> isnscf(x), job.calculations)

cell(job::DFJob)         = cell(structure(job))
calculations(job::DFJob) = job.calculations
isarchived(job::DFJob) = occursin(".archived", job.dir)

"""
    main_job_dir(dir::AbstractString)
    main_job_dir(job::DFJob)

Returns the main directory of the job, also when the job's version is not the one
in the main directory.
"""
main_job_dir(dir::AbstractString) = split(dir, VERSION_DIR_NAME)[1]
main_job_dir(job::DFJob) = main_job_dir(job.dir)

"""
    joinpath(job::DFJob, args...)

`joinpath(job.dir, args...)`.
"""
Base.joinpath(job::DFJob, args...) = joinpath(job.dir, args...)

function Base.pop!(job::DFJob, name::String)
    i = findfirst(x -> x.name == name, job.calculations)
    if i === nothing
        error("No calculation with name $name found.")
    end
    out = job.calculations[i]
    deleteat!(job.calculations, i)
    return out
end

function Utils.searchdir(job::DFJob, str::AbstractString)
    return joinpath.((job,), searchdir(job.dir, str))
end

