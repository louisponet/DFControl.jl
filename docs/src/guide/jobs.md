# Jobs

```@docs
DFJob
```

## Interacting with inputs
```@docs
Base.getindex(::DFJob, ::String)
```
```@setup job_input_access
using DFControl
job = DFJob(joinpath(pathof(DFControl), "..", "..", "docs", "src", "assets", "job"))
```
Example: 
```@repl job_input_access
job["scf"]
job[2]
```

```@docs
Base.push!(::DFJob, ::DFInput)
Base.append!(::DFJob, ::Any...)
Base.pop!(::DFJob)
Base.insert!(::DFJob, ::Int, ::DFInput)
```

## Scheduling, submission and monitoring
```@docs
set_flow!
save(::DFJob)
submit
isrunning
last_running_input
abort
```

## Management

### Directories
```@docs
set_localdir!(::DFJob, ::String)
Base.cp(::DFJob, ::String)
Base.mv(::DFJob, ::String)
Base.filesize
cleanup
```
### Versioning
