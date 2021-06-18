include(joinpath(depsdir, "wannier90flags.jl"))
const WAN_FLAGS = _WAN_FLAGS()
flagtype(::Type{Wannier90}, flag) = haskey(WAN_FLAGS, flag) ? WAN_FLAGS[flag] : Nothing
flagtype(::DFCalculation{Wannier90}, flag) = flagtype(Wannier90, flag)
