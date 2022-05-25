const Point{N,F} = SVector{N,F}
const Point3{F}  = SVector{3,F}
const Vec3{F}    = SVector{3,F}
const Mat3{F}    = SMatrix{3,3,F,9}
const Mat4{F}    = SMatrix{4,4,F,16}
const Vec{N,T}   = SVector{N,T}

function Base.convert(::Type{Point3{T}}, x::Vector{T}) where {T<:AbstractFloat}
    return Point3{T}(x[1], x[2], x[3])
end

abstract type AbstractBand end

"""
Energy band from DFT calculation.
"""
mutable struct Band <: AbstractBand
    k_points_cart  :: Vector{Vec3}
    k_points_cryst :: Vector{Vec3{Float64}}
    eigvals        :: Vector{Float64}
    extra          :: Dict{Symbol,Any}
end
function Band(k_points_cart, k_points_cryst, eigvals)
    return Band(k_points_cart, k_points_cryst, eigvals, Dict{Symbol,Any}())
end

function kpoints(band::Band, kind = :cryst)
    return kind == :cart ? band.k_points_cart : band.k_points_cryst
end
eigvals(band::Band) = band.eigvals
StructTypes.StructType(::Type{<:AbstractBand}) = StructTypes.Mutable()

"""
    bandgap(bands::AbstractVector{Band}, fermi=0.0)

Calculates the bandgap (possibly indirect) around the fermi level.
"""
function bandgap(bands::Union{Iterators.Flatten,AbstractVector{<:AbstractBand}},
                 fermi = 0.0)
    max_valence = -Inf
    min_conduction = Inf
    for b in bands
        max = maximum(eigvals(b) .- fermi)
        min = minimum(eigvals(b) .- fermi)
        if min <= 0.0 <= max
            return 0.0
        end
        if max_valence <= max <= 0.0
            max_valence = max
        end
        if 0.0 <= min <= min_conduction
            min_conduction = min
        end
    end
    return min_conduction - max_valence
end
function bandgap(u_d_bands::Union{NamedTuple,Tuple}, args...)
    return bandgap(Iterators.flatten(u_d_bands), args...)
end

mutable struct TimingData
    name::String
    cpu::Dates.CompoundPeriod
    wall::Dates.CompoundPeriod
    calls::Int
    children::Vector{TimingData}
end
StructTypes.StructType(::Type{TimingData}) = StructTypes.Mutable()
StructTypes.StructType(::Type{Dates.CompoundPeriod}) = StructTypes.Struct()
StructTypes.StructType(::Type{Dates.Millisecond}) = StructTypes.Struct()
StructTypes.StructType(::Type{Dates.Second}) = StructTypes.Struct()
