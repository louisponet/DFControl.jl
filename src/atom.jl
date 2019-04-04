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
open(joinpath(@__DIR__, "..", "assets", "elements.txt"), "r") do f
    while !eof(f)
        line = split(readline(f))
        push!(ELEMENTS, Element(Symbol(line[4]), Meta.parse(line[1]), line[9], Meta.parse(line[10]), (Meta.parse.(line[5:7],)...,)./65535))
    end
end
#undefined element
push!(ELEMENTS, Element(:undef, 0, "undef",0, (0.0, 0.0, 0.0)))
function element(sym::Symbol)
    if tryparse(Int, String(sym)[end:end]) != nothing
        sym = Symbol(String(sym)[1:end-1])
    end
    found = filter(x->x.symbol == sym, ELEMENTS)
    if isempty(found)
       @warn "No element found with symbol '$sym'."
       return ELEMENTS[end]
    end
    return found[1]
end

element(z::Int) = getfirst(x->x.Z == z, ELEMENTS)

abstract type AbstractAtom{T} end

#Definition of the AbstractAtom interface, each AbstractAtom needs to provide a method `atom(...)` that returns a datastructure with the Atom fields.
elsym(atom::AbstractAtom) = element(atom).symbol
"Extracts all the positions of the atoms and puts them in a vector."
positions(atoms::Vector{<:AbstractAtom}, name::Symbol) = [position(x) for x in filter(y -> name(y) == name, atoms)]
getpseudoset(at::AbstractAtom) = getpseudoset(elsym(at), pseudo(at))

"Takes a Vector of atoms and returns a Vector with the atoms having unique symbols."
function Base.unique(atoms::Vector{<:AbstractAtom{T}}) where T <: AbstractFloat
    names    = Symbol[]
    uni = AbstractAtom{T}[]
    for at in atoms
        if !in(name(at), names)
            push!(names, name(at))
            push!(uni, at)
        end
    end
    return uni
end

#We use angstrom everywhere
@with_kw mutable struct Atom{T<:AbstractFloat} <: AbstractAtom{T}
    name        ::Symbol
    element     ::Element
    position    ::Point3{T}
    pseudo      ::String = ""
    projections ::Vector{Projection} = Projection[]
end

Atom(name::Symbol, el::Element, pos::Point3; kwargs...)  = Atom(name=name, element=el, position=pos; kwargs...)
Atom(name::Symbol, el::Symbol, pos::Point3; kwargs...)  = Atom(name=name, element=element(el), position=pos; kwargs...)
Atom(orig_at::Atom, new_pos::Point3) = Atom(name(orig_at), element(orig_at), new_pos, pseudo(orig_at), projections(orig_at))
#Easiest way to implement a new abstractatom is to provide a way to access
#the struct holding `name`, `position`, `element`, `pseudo`, `projection` fields
atom(at::Atom) = at
name(at::AbstractAtom)          = atom(at).name
position(at::AbstractAtom)    = atom(at).position
element(at::AbstractAtom)     = atom(at).element
pseudo(at::AbstractAtom)      = atom(at).pseudo
projections(at::AbstractAtom) = atom(at).projections

setname!(at::AbstractAtom, name::Symbol) = (atom(at).name = name)
setposition!(at::AbstractAtom{T}, position::Point3) where T =
    (atom(at).position = convert(Point3{T}, position))
setpseudo!(at::AbstractAtom, pseudo) = (atom(at).pseudo = pseudo)
setprojections!(at::AbstractAtom, projections::Vector{Projection}) =
    (atom(at).projections = projections)

bondlength(at1::AbstractAtom{T}, at2::AbstractAtom{T}, R=T(0.0)) where T<:AbstractFloat = norm(position(at1) - position(at2) - R)

import Base: ==
==(at1::AbstractAtom, at2::AbstractAtom) =
    name(at1)==name(at2) && norm(position(at1) - position(at2)) < 1e-6
