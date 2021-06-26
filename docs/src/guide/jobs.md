# Jobs
## Contents
```@contents
Pages=["jobs.md"]
```
## Index
```@index
Pages=["jobs.md"]
```
## DFJob
```@docs
DFJob
```

## Interacting with calculations
```@docs
Base.getindex(::DFJob, ::String)
```
```@setup job_calculation_access
using DFControl
job = DFJob(joinpath(pathof(DFControl), "..", "..", "docs", "src", "assets", "job"))
```
Example: 
```@repl job_calculation_access
job["scf"]
job[2]
```

```@docs
Base.push!(::DFJob, ::DFCalculation)
Base.append!(::DFJob, ::Any...)
Base.pop!(::DFJob)
Base.insert!(::DFJob, ::Int, ::DFCalculation)
```

## Scheduling, submission and monitoring
```@docs
set_flow!
save(::DFJob)
submit
isrunning
last_running_calculation
abort
```

## Management

### Directories
```@docs
set_localdir!(::DFJob, ::String)
Base.cp(::DFJob, ::String)
Base.mv(::DFJob, ::String)
Base.filesize
Base.joinpath(::DFJob, ::Any...)
cleanup
```
### Registry
All [`DFJobs`](@ref DFJob) are stored in an internal registry the first time `save(job)` is called. 
This means that finding all previously worked on [`DFJobs`](@ref DFJob) is as straightforward as
calling `DFJob(fuzzy)` where `fuzzy` is a part of the previously saved [`DFJob`](@ref) `local_dir`. 
This will then show a menu in the REPL with the possible choices and one will be loaded upon choosing.

Loading all [`DFJobs`](@ref DFJob) that contain a given `fuzzy` can be done through [`load_jobs(fuzzy)`](@ref load_jobs).

```@docs
registered_jobs
load_jobs
```
### Versioning

As previously mentioned, a rudimentary implementation of a `DFJob` versioning system is implemented. 
Upon calling `save` on a `DFJob`, if there is already a valid `job` script present in `job.local_dir`, 
it is assumed that this was a previous version of the `job` and the script together with all other
files in `job.local_dir` will be copied to a subdirectory of the `.versions` directory bearing the name of the 
respective previous job version. After this, `job.version` will be incremented by `1` signalling the 
new version of the current `job`. 

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
More functionality related to versioning will be added in the future.
