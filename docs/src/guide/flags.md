# Calculation Flags

```@docs
Base.getindex(::DFJob, ::Symbol)
Base.setindex!(::DFJob, ::Any, ::Symbol)
```
```@repl
using DFControl #hide
job = DFJob(joinpath(pathof(DFControl), "..", "..", "docs", "src", "assets", "job")) #hide
job[:ecutwfc]
job[:ecutwfc] = 40
```
