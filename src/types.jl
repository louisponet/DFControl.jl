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
