module Structures
# This module handles all functionality related to Structure and Atom
using ..DFControl
using ..Utils
using LinearAlgebra, StructTypes, Parameters, StaticArrays, JSON3
using spglib_jll
using JLD2
const SPGLIB = spglib_jll.libsymspg

using UnitfulAtomic.Unitful: angstrom, Length, @unit, FreeUnits, unit, 𝐋, FreeUnits,
                             Quantity, ustrip, uconvert
const Ang = angstrom

const ReciprocalType{T,A} = Quantity{T,𝐋^-1,FreeUnits{A,𝐋^-1,nothing}}

StructTypes.StructType(::Type{<:Quantity}) = StructTypes.Struct()
@inline function StaticArrays._inv(::StaticArrays.Size{(3, 3)},
                                   A::SMatrix{3,3,LT}) where {LT<:Union{Length,
                                                                        ReciprocalType}}
    @inbounds x0 = SVector{3}(A[1], A[2], A[3])
    @inbounds x1 = SVector{3}(A[4], A[5], A[6])
    @inbounds x2 = SVector{3}(A[7], A[8], A[9])

    y0 = cross(x1, x2)
    d  = StaticArrays.bilinear_vecdot(x0, y0)
    x0 = x0 / d
    y0 = y0 / d
    y1 = cross(x2, x0)
    y2 = cross(x0, x1)

    @inbounds return SMatrix{3,3}((y0[1], y1[1], y2[1], y0[2], y1[2], y2[2], y0[3], y1[3],
                                   y2[3]))
end

include("projections.jl")
include("atom.jl")
include("structure.jl")

export Projection, Atom, Element, Structure, DFTU
export angstrom, Ang

end
