module Wannier90
    using ..DFControl
    import ..DFControl: Package
    struct Wan90 <: Package end
    include("constants.jl")
    include("input.jl")
    include("fileio.jl")
end
