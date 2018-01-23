

mutable struct Structure{T <: AbstractFloat}
    name ::String
    cell ::Matrix{T}
    atoms::Array{Atom{T}, 1}
    data ::Dict{Symbol, Any}
end

Structure(name, cell::Matrix{T}, atoms::Array{Atom{T}, 1}) where T <: AbstractFloat = Structure(name, cell, atoms, Dict{Symbol, Any}())
Structure(cell::Matrix{T}, atoms::Array{Atom{T}, 1}) where T <: AbstractFloat = Structure("NoName", cell, atoms, Dict{Symbol, Any}())
Structure() = Structure("NoName", eye(3), Atom[], Dict{Symbol, Any}())