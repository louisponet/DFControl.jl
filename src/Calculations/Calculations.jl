module Calculations
    # This module handles all interactions with calculations
    using Parameters, StructTypes, LinearAlgebra
    using ..DFControl
    using ..Utils
    using ..Structures
    
    include("execs.jl")
    include("calculation.jl")
    include("qe.jl")
    include("elk.jl")
    include("wannier.jl")
    include("abinit.jl")

    export Exec, Calculation, InputData
    export Wannier90, QE, Abinit, Elk
end
