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

ismagnetic(at::AbstractAtom) = !iszero(sum(magnetization(at)))

"Takes a Vector of atoms and returns a Vector with the atoms having unique symbols."
function Base.unique(atoms::Vector{A}) where A <: AbstractAtom
    uni = A[]
    for at in atoms
        if findfirst(x -> isequal_species(x, at), uni) === nothing
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
	J::Vector{T} = T[zero(T)]
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
isdefault(x::AbstractVector) = isempty(x) || all(iszero, x)

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
@with_kw_noshow mutable struct Atom{T<:AbstractFloat, LT<:Length{T}} <: AbstractAtom{T, LT}
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
Atom(orig_at::Atom, new_pos_cart::Point3, new_pos_cryst::Point3) = Atom(name(orig_at), element(orig_at), new_pos_cart, new_pos_cryst, pseudo(orig_at), projections(orig_at), magnetization(orig_at), dftu(orig_at))
#Easiest way to implement a new abstractatom is to provide a way to access
#the struct holding `name`, `position_cart`, `element`, `pseudo`, `projection` fields
atom(at::Atom) = at

for interface_function in fieldnames(Atom)
	@eval $interface_function(at::AbstractAtom) = atom(at).$interface_function
	@eval export $interface_function
end

function Base.range(at::AbstractAtom)
	projs = projections(at)
	@assert length(projs) != 0 "At $(name(at)) has no defined projections. Please use `setprojections!` first."
	return projs[1].start : projs[end].last
end

Base.range(v::Vector{AbstractAtom}) = vcat(range.(v)...)

setname!(at::AbstractAtom, name::Symbol) =
	atom(at).name = name

length_unit(at::Atom{T, LT}) where {T, LT} = LT

function setpseudo!(at::AbstractAtom, pseudo::Pseudo; print=true)
	print && @info "Pseudo of atom $(name(at)) set to $pseudo."
	!ispath(path(pseudo)) && @warn "Pseudopath $(path(pseudo)) not found."
	atom(at).pseudo = pseudo
end

function setprojections!(at::AbstractAtom, projections::Vector{Projection}; print=true)
    print && @info "Setting projections for atom $(name(at)) to $projections"
    atom(at).projections = projections
end

bondlength(at1::AbstractAtom{T}, at2::AbstractAtom{T}, R=T(0.0)) where T<:AbstractFloat = norm(position_cart(at1) - position_cart(at2) - R)

==(at1::AbstractAtom{T, LT}, at2::AbstractAtom{T, LT}) where {T, LT} =
    name(at1) == name(at2) && norm(position_cart(at1) - position_cart(at2)) < LT(1e-6)

function isequal_species(at1::AbstractAtom, at2::AbstractAtom)
    for f in fieldnames(typeof(at1))
        if f in (:position_cart, :position_cryst, :projections)
            continue
        end
        if getfield(at1, f) != getfield(at2, f)
            return false
        end
    end
    return true
end

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

function set_magnetization!(at::AbstractAtom, mag; print=true)
	at.magnetization = convert(Vec3, mag)
	print && @info "Magnetization of at $(name(at)) was set to $(magnetization(at))"
end

"""
	distance(at1::AbstractAtom, at2::AbstractAtom)

Calculates the distance between the two atoms.
"""
distance(at1::AbstractAtom, at2::AbstractAtom) = norm(at1.position_cart - at2.position_cart)


"""
	set_position!(at::AbstractAtom, pos::AbstractVector{T}, unit_cell::Mat3) where {T<:Real}

Updates the position of the atom to this. The unit cell is used to make sure both `position_cryst` and `position_cart` are correct.
"""
function set_position!(at::AbstractAtom, pos::AbstractVector{T}, unit_cell::Mat3) where {T<:Real}
	atom(at).position_cryst = Point3{T}(pos...)
	atom(at).position_cart  = unit_cell * atom(at).position_cryst
	return at
end

function set_position!(at::AbstractAtom, pos::AbstractVector{T}, unit_cell::Mat3) where {T<:Length}
	atom(at).position_cart  = Point3{T}(pos...)
	atom(at).position_cryst = unit_cell^-1 * atom(at).position_cart
	return at
end

"""
	scale_bondlength!(at1::AbstractAtom, at2::AbstractAtom, scale::Real, cell::Mat3)

Scales the bondlength between two atoms. The center of mass remains the same.
"""
function scale_bondlength!(at1::AbstractAtom, at2::AbstractAtom, scale::Real, cell::Mat3)
	p1 = position_cryst(at1) 
	p2 = position_cryst(at2) 
	mid = (p1 + p2)/2
	orig_dist = norm(p1 - p2)
	direction = (p1 - mid)/norm(p1-mid)
	new_dist = orig_dist * scale
	new_p1 = mid + new_dist * direction/2
	new_p2 = mid - new_dist * direction/2
	set_position!(at1, new_p1, cell)
	set_position!(at2, new_p2, cell)
end

"""
    polyhedron(at::AbstractAtom, atoms::Vector{<:AbstractAtom}, order::Int)
    polyhedron(at::AbstractAtom, str::AbstractStructure, order::Int)

Returns a polyhedron around the atom, i.e. the `order` closest atoms.
The returned atoms will be ordered according to their distance to the first one.
In the case of a structure rather than a set of atoms, the search will
be performed over all atoms in the structure.
"""
function polyhedron(at::AbstractAtom, atoms::Vector{<:AbstractAtom}, order::Int)
    return sort(atoms, by = x -> distance(x, at))[1:order]
end
