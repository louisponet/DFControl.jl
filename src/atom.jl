include("wannier90/projections.jl")

"""
Represents an element.
"""
struct Element
    symbol::Symbol
    Z::Int64
    name::String
    atomic_weight::Float64
end

"""
Reads all the elements from the file.
"""
const ELEMENTS = Element[]
open(joinpath(@__DIR__, "../assets/elements.txt"), "r") do f
    while !eof(f)
        line = split(readline(f))
        push!(ELEMENTS, Element(Symbol(line[4]), parse(line[1]), line[9], parse(line[10])))
    end
end

function element(sym::Symbol)
    if !isnull(tryparse(Int, String(sym)[end:end]))
        sym = Symbol(String(sym)[1:end-1])
    end
    found = filter(x->x.symbol == sym, ELEMENTS)
    if isempty(found)
        error("No element found with symbol '$sym'.")
    end
    return found[1]
end

function element(z::Int)
    found = filter(x->x.Z == z, ELEMENTS)
    if isempty(found)
        error("No element found with Z '$z'.")
    end
    return found[1]
end

#We use angstrom everywhere
struct Atom{T <: AbstractFloat}
    id      ::Symbol
    element ::Element
    position::Point3D{T}
    data    ::Dict{Symbol, Any}
end

Atom(id::Symbol, element::Element, position::Point3D) = Atom(id, element, position, Dict{Symbol, Any}())
Atom(id::Symbol, element::Symbol, position::Point3D)  = Atom(id, ELEMENTS[element], position, Dict{Symbol, Any}())

positions(atoms::Array{<:Atom, 1}, id::Symbol) = [x.position for x in filter(y -> y.id == id, atoms)]

function unique_atoms(atoms::Array{Atom{T},1 }) where T <: AbstractFloat
    ids    = Symbol[]
    unique = Atom{T}[]
    for at in atoms
        if !in(at.id, ids) 
            push!(ids, at.id)
            push!(unique, at)
        end
    end
    return unique
end

function convert_2atoms(atoms::Array{<:AbstractString, 1}, U=Float64; pseudo_set=:default, pseudo_specifier="")
    out_atoms = Atom{U}[]
    for line in atoms
        atsym, x, y, z = parse.(split(line))
        el = element(atsym)
        pseudo = pseudo_set == :default ? "" : get_default_pseudo(atsym, pseudo_set, pseudo_specifier=pseudo_specifier)
        if pseudo == "" || pseudo == nothing
            push!(out_atoms, Atom{U}(atsym, el, Point3D{U}(x, y, z), Dict{Symbol, Any}()))
        else
            push!(out_atoms, Atom{U}(atsym, el, Point3D{U}(x, y, z), Dict{Symbol, Any}(:pseudo => pseudo)))
        end
    end
    return out_atoms
end

function convert_2atoms(atoms, U=Float64; pseudo_set=:default, pseudo_specifier="")
    out_atoms = Atom{U}[]
    for (atsym, at) in atoms
        el = element(atsym)
        pseudo = pseudo_set == :default ? "" : get_default_pseudo(atsym, pseudo_set, pseudo_specifier=pseudo_specifier)
        if pseudo == "" || pseudo == nothing
            push!(out_atoms, Atom{U}(atsym, el, Point3D{U}(x, y, z), Dict{Symbol, Any}()))
        else
            push!(out_atoms, Atom{U}(atsym, el, Point3D{U}(x, y, z), Dict{Symbol, Any}(:pseudo => pseudo)))
        end
    end
    return out_atoms
end