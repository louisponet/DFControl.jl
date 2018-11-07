module QuantumEspresso

    using ..DFControl
    import ..DFControl: Package
    struct QE <: Package end
    include("constants.jl")
    include("input.jl")
    include("fileio.jl")
    export read_output, read_input

end
