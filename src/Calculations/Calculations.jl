module Calculations
    # This module handles all interactions with calculations
    using ..Utils
    using Parameters
    using StructTypes
    include("execs.jl")
    include("calculation.jl")
    include("qe.jl")

    export Exec, Calculation, InputData
    export Wannier90, QE, Abinit, Elk
end
