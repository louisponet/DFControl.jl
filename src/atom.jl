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
mutable struct Atom{T <: AbstractFloat}
    id          ::Symbol
    element     ::Element
    position    ::Point3D{T}
    pseudo      ::String
    l_soc       ::T
    projections ::Array{Projection, 1}
    mag_moment  ::Array{T, 1}        
    function Atom(id::Symbol, element::Element, position::Point3D{T}, args...) where T <: AbstractFloat
        atom          = new{T}()
        atom.id       = id
        atom.element  = element
        atom.position = position
        names = fieldnames(Atom)
        types = fieldtype.(Atom, names)
        for (name, typ) in zip(names[4:end], types[4:end])
            found = false
            for (field, value) in args
                if field == name
                    setfield!(atom, field, convert(typ, value))
                    found = true
                end
            end
            if !found
                try
                    setfield!(atom, field, zero(typ))
                end
            end
        end
        return atom
    end
end
Atom(id::Symbol, element::Symbol, position::Point3D, args...)  = Atom(id, ELEMENTS[element], position, args...)

positions(atoms::Array{<:Atom, 1}, id::Symbol) = [x.position for x in filter(y -> y.id == id, atoms)]

function unique_atoms(atoms::Array{Atom{T}, 1}) where T <: AbstractFloat
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

function convert_2atoms(atoms::Array{<:AbstractString, 1}, T=Float64)
    out_atoms = Atom{T}[]
    for line in atoms
        atsym, x, y, z = parse.(split(line))
        el = element(atsym)
        push!(out_atoms, Atom(atsym, el, Point3D{U}(x, y, z)))
    end
    return out_atoms
end

function convert_2atoms(atoms, T=Float64)
    out_atoms = Atom{T}[]
    for (atsym, at) in atoms
        el = element(atsym)
        for pos in at
            println(pos)
            push!(out_atoms, Atom(atsym, el, pos))
        end
    end
    return out_atoms
end