const Point{N,F} = SVector{N,F}
const Point3{F} = SVector{3,F}
const Vec3{F} = SVector{3,F}

const Mat3{F} = SMatrix{3,3,F,9}
const Mat4{F} = SMatrix{4,4,F,16}

const SymAnyDict = Dict{Symbol,Any}
const Vec{N,T} = SVector{N,T}

Point{N,T}(x::T) where {N,T} = Point{N,T}(x, x, x)

function Base.convert(::Type{Point3{T}}, x::Vector{T}) where {T<:AbstractFloat}
    return Point3{T}(x[1], x[2], x[3])
end

export Point, Point3, Vec3, Mat3, Mat4
