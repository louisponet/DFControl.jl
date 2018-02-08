
abstract type AbstractStructure{T} end

mutable struct Structure{T <: AbstractFloat} <: AbstractStructure{T}
    name ::AbstractString
    cell ::Mat3{T}
    atoms::Vector{<:AbstractAtom{T}}
    data ::Dict{Symbol, Any}
end

Structure(name, cell::Mat3{T}, atoms::Vector{Atom{T}}) where T <: AbstractFloat = Structure{T}(name, cell, atoms, Dict{Symbol, Any}())
Structure(cell::Matrix{T}, atoms::Vector{Atom{T}}) where T <: AbstractFloat = Structure{T}("NoName", cell, atoms, Dict{Symbol, Any}())
Structure() = Structure("NoName", eye(3), Atom[], Dict{Symbol, Any}())
