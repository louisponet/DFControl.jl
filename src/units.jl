Unitful.register(@__MODULE__)
Base.eltype(::Type{Length{T}}) where {T} = T
using Unitful: angstrom
const Ang = angstrom
@unit e₀ "eₒ" ElementaryCharge 1.602176620898e-19 * u"C" false
@unit kₑ "kₑ" CoulombForceConstant 1 / (4π)u"ϵ0" false
@unit a₀ "a₀" BohrRadius 1u"ħ^2/(1kₑ*me*e₀^2)" false
@unit Eₕ "Eₕ" HartreeEnergy 1u"me*e₀^4*kₑ^2/(1ħ^2)" true
@unit Ry "Ry" RydbergEnergy 0.5Eₕ true

const localunits = Unitful.basefactors
const ReciprocalType{T,A} = Quantity{T,𝐋^-1,FreeUnits{A,𝐋^-1,nothing}}

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
