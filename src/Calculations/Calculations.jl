module Calculations
    # This module handles all interactions with calculations
    using Parameters, StructTypes
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
    
    function Base.hash(data::T, h::UInt) where {T<:Union{<:InputData,Exec}}
        for f in fieldnames(T)
            h = hash(getfield(data, f), h)
        end
        return h
    end

end
