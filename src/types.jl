const SymAnyDict = Dict{Symbol, Any}

Base.convert(::Type{Point3{T}}, x::Vector{T}) where T<:AbstractFloat = Point3{T}(x[1], x[2], x[3])


abstract type Band end

"""
Energy band from DFT calculation.
"""
mutable struct DFBand{T<:AbstractFloat} <: Band
    k_points_cart  ::Vector{Vector{T}}
    k_points_cryst ::Vector{Vector{T}}
    eigvals        ::Vector{T}
    extra          ::Dict{Symbol, Any}
end
DFBand(k_points_cart, k_points_cryst, eigvals) = DFBand(k_points_cart, k_points_cryst, eigvals, Dict{Symbol,Any}())

mutable struct Exec
    exec ::String
    dir  ::String
    flags::Dict{Symbol, Any}
end

Exec(exec::String) = Exec(exec, "~/bin/", SymAnyDict())
Exec(exec::String, dir::String) = Exec(exec, dir, SymAnyDict())
Exec(exec::String, dir::String, flags...) = Exec(exec, dir, SymAnyDict(flags))
