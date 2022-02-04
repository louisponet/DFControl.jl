# Jobs
## Contents
```@contents
Pages=["jobs.md"]
```
## Index
```@index
Pages=["jobs.md"]
```
## Job
```@docs
Job
load(::Server, ::Job)
```

## Interacting with calculations
```@docs
Base.getindex(::Job, ::String)
```
```@setup job_calculation_access
using DFControl
job = load(Job(joinpath(@__DIR__, "..", "assets", "job")))
```
Example: 
```@repl job_calculation_access
job["scf"]
job[2]
```

```@docs
Base.push!(::Job, ::Calculation)
Base.append!(::Job, ::Any...)
Base.pop!(::Job)
Base.insert!(::Job, ::Int, ::Calculation)
```

## Scheduling, submission and monitoring
```@docs
set_flow!
save(::Job)
submit
isrunning
abort
```

## Directories
```@docs
joinpath(::Job, ::Any...)
abspath(::Job)
cleanup(::Job)
```
## Registry
All [`Jobs`](@ref Job) are stored in an internal registry the first time `save(job)` is called. 
This means that finding all previously worked on [`Jobs`](@ref Job) is as straightforward as
calling `load(server, Job(fuzzy))` where `fuzzy` is a part of the previously saved [`Job`](@ref) `dir`. 
This will then return a list of [`Jobs`](@ref Job) with similar directories. 

## Versioning

As previously mentioned, a rudimentary implementation of a [`Job`](@ref) versioning system is implemented. 
Upon calling `save` on a [`Job`](@ref), if there is already a valid `job` script present in `job.dir`, 
it is assumed that this was a previous version of the `job` and the script together with all other
files in `job.local_dir` will be copied to a subdirectory of the `.versions` directory bearing the name of the 
respective previous job version. After this, `job.version` will be incremented by `1` signalling the 
new version of the current [`Job`](@ref). 

The virtue of this system is that it is possible to roll back to a previous version after possibly making
breaking changes, or to pull out previous results after further experimentation was performed.

!!! note
    If `job.copy_temp_folders=true` all possible intermediate files inside the temporary calculation directory 
    (i.e. "job_dir/outputs") will be copied every time the `job` is saved. These can be quite large and 
    can quickly create very large `job` directories. Handle with care!

```@docs
versions
last_version
switch_version!
rm_version!
```

## Archiving
After a `Job` is completed, or an interesting result is achieved, it makes sense to store it for future reference. 
This can be achieved through the [`archive`](@ref) function. This will take the current job, and copy it to a subdirectory (specified by the second argument to [`archive`](@ref)) of the `jobs/archived` directory inside the `DFControl` config directory. The third argument is a description of this job's result.

!!! note
    In order to not cause huge file transfers, all the temporary directories will first be removed before archiving.

Example:
```julia
archive(job, "test_archived_job", "This is a test archived job")
```

To query previously archived jobs one can use `load(Server("localhost"), Job("archived"))`.

```@docs
archive
```

## Output 
```@docs
outputdata
readfermi
readbands
bandgap
```

## [Environments](@id environments_header)
Environments specify the skeleton of the job script, i.e. which environment variables need to be set, which scheduler flags, etc.

```@docs
Environment
```
