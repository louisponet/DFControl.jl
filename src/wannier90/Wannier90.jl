module Wannier90
    using ..DFControl
    using ..DFControl.GeometryTypes
    using ..DFControl.LinearAlgebra
    import ..DFControl: Package, InputData, DFInput, Projection, Atom, Structure
    struct Wan90 <: Package end
    include("constants.jl")
    include("input.jl")
    include("fileio.jl")
end
