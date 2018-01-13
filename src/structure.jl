#We use angstrom everywhere
struct Atom{T <: AbstractFloat}
    id::Symbol
    element::Element
    position::Point3D{T}
end

mutable struct Structure{T <: AbstractFloat}
    name::String
    cell::Matrix{T}
    atoms::Array{Atom{T}, 1}
end

