module QuantumEspresso

    using ..DFControl
    using ..DFControl.GeometryTypes
    using ..DFControl.LinearAlgebra
    import ..DFControl: Package, InputData, DFInput, Atom, Structure, DFJob, DFBand

    struct QE <: Package end
    include("constants.jl")
    include("input.jl")
    include("fileio.jl")
    export read_output, read_input

end
