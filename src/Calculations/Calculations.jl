module Calculations
    # This module handles all interactions with calculations
    using ..Utils
    using Parameters
    using StructTypes

    include("qe.jl")
    include("wannier.jl")
    
end
