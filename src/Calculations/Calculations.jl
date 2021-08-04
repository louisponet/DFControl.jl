module Calculations
    # This module handles all interactions with calculations
    using ..Utils
    using ..DFControl
    using Parameters
    using StructTypes
    include("execs.jl")
    include("calculation.jl")
    include("qe.jl")
    include("elk.jl")
    include("wannier.jl")
    include("abinit.jl")
    

    export Exec, Calculation, InputData
    export Wannier90, QE, Abinit, Elk
end
