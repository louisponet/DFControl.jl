include("wannier90/projections.jl")

import Base: ==
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

@with_kw struct DFTU{T}
	l::Int = -1
	U::T   = zero(T)
	J0::T  = zero(T)
	#QE params
	α::T   = zero(T)
	β::T   = zero(T)
	J::Vector{T} = T[]
end

function ==(x::DFTU, y::DFTU)
	fnames = fieldnames(DFTU)
	for fn in fnames
		if getfield(x, fn) != getfield(y, fn)
			return false
		end
	end
	return true
end

isdefault(x::DFTU{T}) where {T} =
	x == DFTU{T}() 

isdefault(x::Any) = isempty(x)

Base.string(::Type{Elk}, dftu::DFTU) = "$(dftu.l) $(dftu.U) $(dftu.J0)"

# TODO Multiple l per atom in Elk??
#We use angstrom everywhere
@with_kw mutable struct Atom{T<:AbstractFloat} <: AbstractAtom{T}
    name          ::Symbol
    element       ::Element
    position      ::Point3{T}
    pseudo        ::String = ""
    projections   ::Vector{Projection} = Projection[]
    magnetization ::Vec3{T} = zero(Vec3{T})
    dftu          ::DFTU{T} = DFTU{T}()
end

Atom(name::Symbol, el::Element, pos::Point3{T}; kwargs...) where {T} = Atom{T}(name=name, element=el, position=pos; kwargs...)
Atom(name::Symbol, el::Symbol, pos::Point3{T}; kwargs...) where {T} = Atom{T}(name=name, element=element(el), position=pos; kwargs...)
Atom(orig_at::Atom, new_pos::Point3) = Atom(name(orig_at), element(orig_at), new_pos, pseudo(orig_at), projections(orig_at), magnetization(orig_at), dftu(orig_at))
#Easiest way to implement a new abstractatom is to provide a way to access
#the struct holding `name`, `position`, `element`, `pseudo`, `projection` fields
atom(at::Atom) = at

for interface_function in (:name, :position, :element, :pseudo, :projections, :magnetization, :dftu)
	@eval $interface_function(at::AbstractAtom) = atom(at).$interface_function
end


function Base.range(at::AbstractAtom)
	projs = projections(at)
	@assert length(projs) != 0 "At $(name(at)) has no defined projections. Please use `setprojections!` first."
	return projs[1].start : projs[end].last
end

setname!(at::AbstractAtom, name::Symbol) =
	atom(at).name = name

setposition!(at::AbstractAtom{T}, position::Point3) where T =
    atom(at).position = convert(Point3{T}, position)

function setpseudo!(at::AbstractAtom, pseudo; print=true)
	print && @info "Pseudo of atom $(name(at)) set to $pseudo."
	atom(at).pseudo = pseudo
end

setprojections!(at::AbstractAtom, projections::Vector{Projection}) =
    atom(at).projections = projections

bondlength(at1::AbstractAtom{T}, at2::AbstractAtom{T}, R=T(0.0)) where T<:AbstractFloat = norm(position(at1) - position(at2) - R)

==(at1::AbstractAtom, at2::AbstractAtom) =
    name(at1) == name(at2) && norm(position(at1) - position(at2)) < 1e-6

#TODO fix documentation
for hub_param in (:U, :J0, :α, :β)
	f = Symbol("set_Hubbard_$(hub_param)!")
	str = "$(hub_param)"
	@eval begin
		"""
			$($(f))(at::AbstractAtom, v::AbstractFloat)

		Set the Hubbard $($(str)) parameter for the specified atom.

		Example:
			`$($(f))(at, 2.1)'
		"""
		function $f(at::AbstractAtom{T}, v::AbstractFloat) where {T}
			dftu(at).$(hub_param) = convert(T, v)
		end
	end
end

"""
	set_Hubbard_J!(at::AbstractAtom, v::Vector{AbstractFloat})

Set the Hubbard J parameter for the specified atom.

Example:
	`set_Hubbard_J(at, [2.1])'
"""
set_Hubbard_J!(at::AbstractAtom{T}, v::AbstractFloat) where {T} =
	dftu(at).J = convert.(T, v)
