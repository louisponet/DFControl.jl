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
        push!(ELEMENTS, Element(Symbol(line[4]), parse(line[1]), line[9], parse(line[10]))
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

