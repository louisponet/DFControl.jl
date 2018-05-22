include("wannier90/projections.jl")

"""
Represents an element.
"""
struct Element
    symbol        ::Symbol
    Z             ::Int64
    name          ::String
    atomic_weight ::Float64
    color         ::NTuple{3, Float64}
end

"""
Reads all the elements from the file.
"""
const ELEMENTS = Element[]
open(joinpath(@__DIR__, "../assets/elements.txt"), "r") do f
    while !eof(f)
        line = split(readline(f))
        push!(ELEMENTS, Element(Symbol(line[4]), parse(line[1]), line[9], parse(line[10]), (parse.(line[5:7])...)./65535))
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

abstract type AbstractAtom{T} end

#Definition of the AbstractAtom interface, each AbstractAtom needs to provide:
#   position(atom)
#   element(atom)
#   id(atom)
#   pseudo(atom)
#   setpseudo!(atom, pseudo)
#   setprojections!(atom, projections)
elsym(atom::AbstractAtom) = element(atom).symbol
"Extracts all the positions of the atoms and puts them in a vector."
positions(atoms::Vector{<:AbstractAtom}, id::Symbol) = [position(x) for x in filter(y -> id(y) == id, atoms)]
getpseudoset(at::AbstractAtom) = getpseudoset(elsym(at), pseudo(at))

"Takes a Vector of atoms and returns a Vector with the atoms having unique symbols."
function unique_atoms(atoms::Vector{<:AbstractAtom{T}}) where T <: AbstractFloat
    ids    = Symbol[]
    unique = AbstractAtom{T}[]
    for at in atoms
        if !in(id(at), ids)
            push!(ids, id(at))
            push!(unique, at)
        end
    end
    return unique
end

#We use angstrom everywhere
mutable struct Atom{T<:AbstractFloat} <: AbstractAtom{T}
    id          ::Symbol
    element     ::Element
    position    ::Point3{T}
    pseudo      ::String
    projections ::Vector{Projection}
    function Atom(id::Symbol, element::Element, position::Point3{T}, args...) where T <: AbstractFloat
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
Atom(id::Symbol, el::Symbol, position::Point3, args...)  = Atom(id, element(el), position, args...)

position(atom::Atom)    = atom.position
element(atom::Atom)     = atom.element
id(atom::Atom)          = atom.id
pseudo(atom::Atom)      = atom.pseudo
projections(atom::Atom) = atom.projections

function setpseudo!(atom::Atom, pseudo)
    atom.pseudo = pseudo
end

function setprojections!(atom::Atom, projections)
    atom.projections = projections
end

function convert_2atoms(atoms::Vector{<:AbstractString}, T=Float64)
    out_atoms = Atom{T}[]
    for line in atoms
        atsym, x, y, z = parse.(split(line))
        el = element(atsym)
        push!(out_atoms, Atom(atsym, el, Point3{U}(x, y, z)))
    end
    return out_atoms
end

function convert_2atoms(atoms, T=Float64)
    out_atoms = Atom{T}[]
    for (atsym, at) in atoms
        el = element(atsym)
        for pos in at
            push!(out_atoms, Atom(atsym, el, pos))
        end
    end
    return out_atoms
end


"Returns the atom to which a certain orbital index belongs."
function orbital2atom(oid, atoms)
    for at in atoms
        for proj in projections(at)
            if proj.start <= oid <= proj.last
                return at
            end
        end
    end
end

bondlength(at1::AbstractAtom{T}, at2::AbstractAtom{T}, R=T(0.0)) where T<:AbstractFloat = norm(position(at1) - position(at2) - R)
