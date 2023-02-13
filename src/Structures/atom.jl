"""
    DFTU(;l ::Int = -1,
          U ::T   = zero(T),
          J0::T  = zero(T),
          α ::T   = zero(T),
          β ::T   = zero(T),
          J ::Vector{T} = T[zero(T)])

DFT+U parameters for a given [`Atom`](@ref).
"""
Base.@kwdef mutable struct DFTU
    l::Int      = -1
    U::Float64  = 0.0
    J0::Float64 = 0.0
    #QE params
    α::Float64 = 0.0
    β::Float64 = 0.0
    J::Vector{Float64} = [0.0]
end
function DFTU(dict::JSON3.Object)
    return DFTU(;dict...)
end

function Base.:(==)(x::DFTU, y::DFTU)
    fnames = fieldnames(DFTU)
    for fn in fnames
        if getfield(x, fn) != getfield(y, fn)
            return false
        end
    end
    return true
end

StructTypes.StructType(::Type{DFTU}) = StructTypes.Struct()

"""
    Element(symbol::Symbol, Z::Int, name::String, atomic_weight::Float64, color::NTuple{3, Float64})
    
Represents an element. Most conveniently used trough the function [`element`](@ref).
"""
struct Element
    symbol        :: Symbol
    Z             :: Int64
    name          :: String
    atomic_weight :: Float64
    color         :: NTuple{3,Float64}
end
Element() = Element(:nothing, -1, "", -1, (0.0, 0.0, 0.0))
StructTypes.StructType(::Type{Element}) = StructTypes.Struct()
function Element(dict::JSON3.Object)
    return Element(Symbol(dict[:symbol]), dict[:Z], dict[:name], dict[:atomic_weight], ([v for v in dict[:color]]...,))
end

"""
Reads all the elements from the file.
"""
const ELEMENTS = Element[]

open(joinpath(@__DIR__, "..", "..", "assets", "elements.txt"), "r") do f
    while !eof(f)
        line = split(readline(f))
        push!(ELEMENTS,
              Element(Symbol(line[4]), Meta.parse(line[1]), line[9], Meta.parse(line[10]),
                      (Meta.parse.(line[5:7],)...,) ./ 65535))
    end
end
#undefined element
push!(ELEMENTS, Element(:undef, 0, "undef", 0, (0.0, 0.0, 0.0)))

"""
    element(sym::Symbol)

Returns the predefined [`Element`](@ref) with symbol `sym`,
i.e. `element(:Si)` will return the pregenerated Silicon [`Element`](@ref).
"""
function element(sym::Symbol)
    # This we do to clean up atom names such as :Fe1
    if tryparse(Int, String(sym)[end:end]) !== nothing
        sym = Symbol(String(sym)[1:end-1])
    end
    found = filter(x -> x.symbol == sym, ELEMENTS)
    if isempty(found)
        return ELEMENTS[end]
    end
    return found[1]
end

element(z::Int) = getfirst(x -> x.Z == z, ELEMENTS)

Base.@kwdef mutable struct Pseudo
    server::String = ""
    path::String = ""
    pseudo::String = ""
end
Pseudo(d::Dict{Symbol, Any}) = Pseudo(; d...)
function Pseudo(d::Dict{String, Any})
    td = Dict{Symbol, Any}()
    for (k, v) in d
        td[Symbol(k)] = v
    end
    Pseudo(td)
end

Base.convert(::Type{Pseudo}, d::Dict) = Pseudo(d)
StructTypes.StructType(::Type{Pseudo}) = StructTypes.Struct()
Base.write(f::AbstractString, p::Pseudo, args...) = write(f, p.pseudo, args...)
Base.write(f::IO, p::Pseudo) = write(f, p.pseudo)
function Base.:(==)(p1::Pseudo, p2::Pseudo)
    if !(isempty(p1.path) && isempty(p2.path) && isempty(p1.server) && isempty(p2.server))
        return p1.path == p2.path && p1.server==p2.server
    else
        return p1.pseudo == p2.pseudo
    end
end
        

# TODO Multiple l per atom in Elk??
#We use angstrom everywhere
"""
    Atom(name::Symbol, element::Element, position_cart::Point3{Length}, position_cryst::Point3;
         pseudo::String = "",
         projections::Vector{Projection} = Projection[],
         magnetization::Vec3 = Vec3(0.0, 0.0, 0.0),
         dftu::DFTU = DFTU())
         
Representation of an `atom`.
    
The `name` of the `atom` is used as an identifier for the `atom` type, in the sense that atoms with the same `pseudo`, `projections`, `magnetization` and [`dftu`](@ref DFTU) attributes should belong to the same type. This also means that during sanity checks atoms that are not of the same type will be given different names. This is done in this way because it often makes sense to change these parameters on all atoms of the same kind at the same time, but still allow the flexibility to change them for individual atoms as well.

`position_cart` should have a valid `Unitful.Length` type such as `Ang`.

See documentation for [`Element`](@ref) for further information on this attribute.
"""
@with_kw_noshow mutable struct Atom
    name::Symbol
    position_cart::Point3{typeof(0.0Ang)}
    position_cryst::Point3{Float64}
    element::Element = element(name)
    pseudo::Pseudo = Pseudo("", "", "")
    projections::Vector{Projection} = Projection[]
    magnetization::Vec3{Float64} = zero(Vec3{Float64})
    dftu::DFTU = DFTU()
end
StructTypes.StructType(::Type{Atom}) = StructTypes.Struct()
function Atom(dict::JSON3.Object)
    return Atom(Symbol(dict[:name]), 1.0Ang .* Point3([t[:val] for t in dict[:position_cart]]),
                Point3([t for t in dict[:position_cryst]]), Element(dict[:element]), dict[:pseudo],
                [Projection(x) for x in dict[:projections]], Vec3([t for t in dict[:magnetization]]),
                DFTU(dict[:dftu]))
end


"Takes a Vector of atoms and returns a Vector with the atoms having unique symbols."
function Base.unique(atoms::Vector{Atom})
    uni = Atom[]
    for at in atoms
        if findfirst(x -> isequal_species(x, at), uni) === nothing
            push!(uni, at)
        end
    end
    return uni
end

function isequal_species(at1::Atom, at2::Atom)
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

function Base.range(at::Atom)
    projs = at.projections
    @assert length(projs) != 0 "At $(at.name) has no defined projections. Please use `setprojections!` first."
    return projs[1].start:projs[end].last
end

Base.range(v::Vector{Atom}) = vcat(range.(v)...)

function projections_string(at::Atom)
    return "$(at.name): " * join(map(x -> x.orbital.name, at.projections), ";")
end

function Base.:(==)(at1::Atom, at2::Atom)
    return at1.name == at2.name && norm(at1.position_cart - at2.position_cart) < 1e-6Ang && at1.dftu == at2.dftu
end

"""
	set_position!(at::Atom, pos::AbstractVector{T}, unit_cell::Mat3) where {T<:Real}

Updates the position of the atom to this. The unit cell is used to make sure both `position_cryst` and `position_cart` are correct.
"""
function set_position!(at::Atom, pos::AbstractVector{T}, unit_cell::Mat3) where {T<:Real}
    at.position_cryst = Point3{T}(pos...)
    at.position_cart  = unit_cell * at.position_cryst
    return at
end

function set_position!(at::Atom, pos::AbstractVector{T}, unit_cell::Mat3) where {T<:Length}
    at.position_cart  = Point3{T}(pos...)
    at.position_cryst = unit_cell^-1 * at.position_cart
    return at
end

"""
	scale_bondlength!(at1::Atom, at2::Atom, scale::Real, cell::Mat3)

Scales the bondlength between two atoms. The center of mass remains the same.
"""
function scale_bondlength!(at1::Atom, at2::Atom, scale::Real, cell::Mat3)
    p1 = at.position_cryst
    p2 = at.position_cryst
    mid = (p1 + p2) / 2
    orig_dist = norm(p1 - p2)
    direction = (p1 - mid) / norm(p1 - mid)
    new_dist = orig_dist * scale
    new_p1 = mid + new_dist * direction / 2
    new_p2 = mid - new_dist * direction / 2
    set_position!(at1, new_p1, cell)
    return set_position!(at2, new_p2, cell)
end

Base.getindex(A::Matrix, a1::T , a2::T) where {T<:Union{Atom, Projection}} =
    getindex(A, range(a1), range(a2))

Base.getindex(A::Matrix, a::Atom) =
    getindex(A, a, a)

Base.getindex(A::Vector, a::Atom) =
    getindex(A, range(a))

Base.view(A::Matrix, a1::T, a2::T) where {T<:Union{Atom, Projection}} =
    view(A, range(a1), range(a2))

Base.view(A::Matrix, a::Union{Atom, Projection}) =
    view(A, range(a), range(a))

Base.view(A::Vector, a::Union{Atom, Projection}) =
    view(A, range(a))

