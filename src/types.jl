const KeyValIterable = Union{Array{Pair{Symbol, Any}, 1}, Dict{Symbol, Any}}
"Point in 3D space in cartesian coordinates with specified float type"
struct Point3D{T<:AbstractFloat}
    x::T
    y::T
    z::T
end
Point3D()                                   = Point3D(0.0)
Point3D(x::T) where T<:AbstractFloat        = Point3D{T}(x, x, x)
Point3D(::Type{T},x) where T<:AbstractFloat = Point3D{T}(x, x, x)
Point3D(x::Array{<:AbstractFloat, 1})        = Point3D(x[1], x[2], x[3])
convert(::Type{Point3D}, x::Array{<:AbstractFloat, 1}) = Point3D(x[1], x[2], x[3])
Point3D{T}() where T<:AbstractFloat         = Point3D{T}(0)

import Base: +, -, *, /, convert, promote_rule, show, zero, norm
+(x::Point3D, y::Point3D) = Point3D(x.x + y.x, x.y + y.y, x.z + y.z)
-(x::Point3D, y::Point3D) = Point3D(x.x - y.x, x.y - y.y, x.z - y.z)
*(x::Point3D, y::Point3D) = Point3D(x.x * y.x, x.y * y.y, x.z * y.z)
*(x::Point3D, y::Number)  = Point3D(x.x * y, x.y * y, x.z * y)
*(y::Number, x::Point3D)  = Point3D(x.x * y, x.y * y, x.z * y)
*(x::Matrix, y::Point3D)  = Point3D(x * Array(y))
/(a::Point3D, b::Point3D) = Point3D(a.x / b.x, a.y / b.y, a.z / b.z)
/(a::Point3D, b::Number)  = Point3D(a.x / b, a.y / b, a.z / b)
@inline norm(a::Point3D)  = sqrt(a.x^2 + a.y^2 + a.z^2)
zero(::Type{Point3D{T}}) where T = Point3D(zero(T), zero(T), zero(T))
zero(::Type{Point3D})     = Point3D(zero(Float64), zero(Float64), zero(Float64))

convert(::Type{Point3D}, x::Point3D) = x
convert(::Type{Array}, x::Point3D)   = [x.x, x.y, x.z]
convert(::Type{Point3D}, x::T) where T<:AbstractFloat             = Point3D{T}(x, x, x)
convert(::Type{Point3D{T}}, x::Real) where T<:AbstractFloat       = Point3D{T}(x, x, x)
convert(::Type{Point3D{T}}, x::Point3D) where T<:AbstractFloat    = Point3D{T}(x.x, x.y, x.z)
convert(::Type{Point3D{T}}, x::Array{T,1}) where T<:AbstractFloat = Point3D{T}(x[1], x[2], x[3])
promote_rule(::Type{Point3D{S}}, ::Type{T}) where {S<:AbstractFloat,T<:Real}                   = Point3D{promote_type(S,T)}
promote_rule(::Type{Point3D{S}}, ::Type{Point3D{T}}) where {S<:AbstractFloat,T<:AbstractFloat} = Point3D{promote_type(S,T)}

show(io::IO, x::Point3D)=print(io, "x = $(x.x), y = $(x.y), z = $(x.z)")
Base.write(f::IO, x::Point3D) = write(f, "$(x.x) $(x.y) $(x.z)")

include("atom.jl")
include("structure.jl")

abstract type Band end

"""
Energy band from DFT calculation.
"""
mutable struct DFBand{T<:AbstractFloat} <: Band
    k_points_cart::Array{Array{T,1},1}
    k_points_cryst::Array{Array{T,1},1}
    eigvals::Array{T,1}
    extra::Dict{Symbol,Any}
end
DFBand(k_points_cart::Array{Array{T,1},1}, k_points_cryst::Array{Array{T,1},1}, eigvals::Array{T,1}) where T <: AbstractFloat = DFBand{T}(k_points_cart, k_points_cryst, eigvals, Dict{Symbol,Any}())


function Base.display(band::DFBand{T}) where T <: AbstractFloat
    string = """DFBand{$T}:
    k_points of length $(length(band.k_points_cryst)):
    cart:    $(band.k_points_cart[1]) -> $(band.k_points_cart[end])
    cryst:   $(band.k_points_cryst[1]) -> $(band.k_points_cryst[end])
    eigvals: $(band.eigvals[1]) -> $(band.eigvals[end])
    extra:   $(band.extra)
    """
    dfprintln(string)
end

function Base.display(bands::Array{<:DFBand})
    map(display,bands)
end