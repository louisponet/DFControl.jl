
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

"""
Returns all the atoms inside the structure with the specified symbol
"""
function get_atoms(str::AbstractStructure, atsym::Symbol)
    out = AbstractAtom[]
    for at in str.atoms
        at.id == atsym && push!(out, at)
    end
    return out
end

"""
Changes the projections of the specified atoms.
"""
function change_projections!(str::Structure, projections...)
    projdict = Dict(projections)
    for at in unique_atoms(str.atoms)
        if !haskey(projdict, at.id)
            projdict[at.id] = [proj.orb for proj in at.projections]
        end
    end
    empty_projections!(str)
    add_projections(projdict, str.atoms)
end

function empty_projections!(str::Structure)
    for at in str.atoms
        empty!(at.projections)
    end
end
