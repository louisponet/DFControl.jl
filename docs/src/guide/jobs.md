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
```

## Interacting with calculations
```@docs
Base.getindex(::Job, ::String)
```
```@setup job_calculation_access
using DFControl
job = Job(joinpath(@__DIR__, "..", "assets", "job"))
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
Jobs.set_flow!
save(::Job)
submit
isrunning
abort
```


## Directories
```@docs
Base.cp(::Job, ::String)
Base.mv(::Job, ::String)
Base.filesize
Base.joinpath(::Job, ::Any...)
```
## Registry
All [`Jobs`](@ref Job) are stored in an internal registry the first time `save(job)` is called. 
This means that finding all previously worked on [`Jobs`](@ref Job) is as straightforward as
calling `Job(fuzzy)` where `fuzzy` is a part of the previously saved [`Job`](@ref) `dir`. 
This will then show a menu in the REPL with the possible choices and one will be loaded upon choosing.

```@docs
registered_jobs
```
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
This can be achieved through the [`Jobs.archive`](@ref) function. This will take the current job, and copy it to a subdirectory (specified by the second argument to [`Jobs.archive`](@ref)) of the `.archived` directory inside the `DFControl` config directory. The third argument is a description of this job's result. If a previous version of the job should be saved, this can be done by specifying the `version` keyword. 

Example:
```julia
archive(job, "test_archived_job", "This is a test archived job"; version = 1, present = j -> @show j.structure)
```

To query previously archived jobs one can use the [`Jobs.archived_jobs`](@ref) function.

```@docs
Jobs.archive
Jobs.archived_jobs
```

## Output 
```@docs
outputdata
readfermi
readbands
bandgap
```

