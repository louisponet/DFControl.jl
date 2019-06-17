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

abstract type AbstractAtom{T, LT<:Length{T}} end

#Definition of the AbstractAtom interface, each AbstractAtom needs to provide a method `atom(...)` that returns a datastructure with the Atom fields.
elsym(atom::AbstractAtom) = element(atom).symbol
"Extracts all the positions of the atoms and puts them in a vector."
positions(atoms::Vector{<:AbstractAtom}, name::Symbol) = [position_cart(x) for x in filter(y -> name(y) == name, atoms)]
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

@with_kw mutable struct DFTU{T}
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

mutable struct Pseudo
	name::String
	dir ::String
	Pseudo() =
		new("", "")
	Pseudo(name::AbstractString, dir::AbstractString) =
		new(name, abspath(dir))
end

Base.isempty(p::Pseudo) =
	isempty(p.name) && isempty(p.dir)

==(p1::Pseudo, p2::Pseudo) =
	p1.name == p2.name && p1.dir == p2.dir

path(p::Pseudo) = joinpath(p.dir, p.name)
# TODO Multiple l per atom in Elk??
#We use angstrom everywhere
@with_kw mutable struct Atom{T<:AbstractFloat, LT<:Length{T}} <: AbstractAtom{T, LT}
    name          ::Symbol
    element       ::Element
    position_cart ::Point3{LT}
    position_cryst::Point3{T}=zero(Point3{T})
    pseudo        ::Pseudo=Pseudo()
    projections   ::Vector{Projection} = Projection[]
    magnetization ::Vec3{T} = zero(Vec3{T})
    dftu          ::DFTU{T} = DFTU{T}()
end

Atom(name::Symbol, el::Element, pos_cart::Point3{LT}, pos_cryst::Point3{T}; kwargs...) where {T, LT<:Length{T}} =
	Atom{T, LT}(name=name, element=el, position_cart=pos_cart, position_cryst=pos_cryst; kwargs...)
Atom(name::Symbol, el::Symbol, args...; kwargs...) =
	Atom(name, element(el), args...; kwargs...)

#TODO this is a little iffy
Atom(orig_at::Atom, new_pos_cart::Point3) = Atom(name(orig_at), element(orig_at), new_pos_cart, position_cryst(orig_at), pseudo(orig_at), projections(orig_at), magnetization(orig_at), dftu(orig_at))
#Easiest way to implement a new abstractatom is to provide a way to access
#the struct holding `name`, `position_cart`, `element`, `pseudo`, `projection` fields
atom(at::Atom) = at

for interface_function in (:name, :position_cart, :position_cryst, :element, :pseudo, :projections, :magnetization, :dftu)
	@eval $interface_function(at::AbstractAtom) = atom(at).$interface_function
end


function Base.range(at::AbstractAtom)
	projs = projections(at)
	@assert length(projs) != 0 "At $(name(at)) has no defined projections. Please use `setprojections!` first."
	return projs[1].start : projs[end].last
end

setname!(at::AbstractAtom, name::Symbol) =
	atom(at).name = name

length_unit(at::Atom{T, LT}) where {T, LT} = LT

setposition_cart!(at::AbstractAtom{T}, position_cart::Point3) where T =
    atom(at).position_cart = length_unit(at).(convert(Point3{T}, position_cart))

function setpseudo!(at::AbstractAtom, pseudo::Pseudo; print=true)
	print && @info "Pseudo of atom $(name(at)) set to $pseudo."
	!ispath(path(pseudo)) && @warn "Pseudopath $(path(pseudo)) not found."
	atom(at).pseudo = pseudo
end

setprojections!(at::AbstractAtom, projections::Vector{Projection}) =
    atom(at).projections = projections

bondlength(at1::AbstractAtom{T}, at2::AbstractAtom{T}, R=T(0.0)) where T<:AbstractFloat = norm(position_cart(at1) - position_cart(at2) - R)

==(at1::AbstractAtom{T, LT}, at2::AbstractAtom{T, LT}) where {T, LT} =
    name(at1) == name(at2) && norm(position_cart(at1) - position_cart(at2)) < LT(1e-6)

#TODO fix documentation
for hub_param in (:U, :J0, :α, :β)
	f = Symbol("set_Hubbard_$(hub_param)!")
	str = "$(hub_param)"
	@eval begin
		"""
			$($(f))(at::AbstractAtom, v::AbstractFloat; print=true)

		Set the Hubbard $($(str)) parameter for the specified atom.

		Example:
			`$($(f))(at, 2.1)'
		"""
		function $f(at::AbstractAtom{T}, v::AbstractFloat; print=true) where {T}
			dftu(at).$(hub_param) = convert(T, v)
			print && @info "Hubbard $($(str)) of atom $(at.name) set to $v"
		end
	end
end

"""
	set_Hubbard_J!(at::AbstractAtom, v::Vector{<:AbstractFloat}; print=true)

Set the Hubbard J parameter for the specified atom.

Example:
	`set_Hubbard_J(at, [2.1])'
"""
function set_Hubbard_J!(at::AbstractAtom{T}, v::Vector{<:AbstractFloat}; print=true) where {T}
	dftu(at).J = convert.(T, v)
	print && @info "Hubbard J of atom $(at.name) set to $v"
end

function position_string(::Type{QE}, at::AbstractAtom; relative=true)
	pos = relative ? round.(position_cryst(at), digits=5) : round.(ustrip.(uconvert.(Ang, position_cart(at))), digits=5)
	return "$(name(at))  $(pos[1]) $(pos[2]) $(pos[3])\n"
end
