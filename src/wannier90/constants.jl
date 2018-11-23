
include(joinpath(depsdir, "wannier90flags.jl"))
const WANFLAGS = _WANFLAGS()
flagtype(::Type{Wannier90}, flag) = haskey(WANFLAGS, flag) ? WANFLAGS[flag] : Nothing
flagtype(::DFInput{Wannier90}, flag) = flagtype(Wannier90, flag)
