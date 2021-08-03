include("projections.jl")

import Base: ==
"""
Reads all the elements from the file.
"""
const ELEMENTS = Element[]
open(joinpath(@__DIR__, "..", "assets", "elements.txt"), "r") do f
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
    if tryparse(Int, String(sym)[end:end]) != nothing
        sym = Symbol(String(sym)[1:end-1])
    end
    found = filter(x -> x.symbol == sym, ELEMENTS)
    if isempty(found)
        @warn "No element found with symbol '$sym'."
        return ELEMENTS[end]
    end
    return found[1]
end

element(z::Int) = getfirst(x -> x.Z == z, ELEMENTS)

#Definition of the AbstractAtom interface, each AbstractAtom needs to provide a method `atom(...)` that returns a datastructure with the Atom fields.
elsym(atom::AbstractAtom) = element(atom).symbol
"Extracts all the positions of the atoms and puts them in a vector."
function positions(atoms::Vector{<:AbstractAtom}, name::Symbol)
    return [position_cart(x) for x in filter(y -> name(y) == name, atoms)]
end
getpseudoset(at::AbstractAtom) = getpseudoset(elsym(at), pseudo(at))

ismagnetic(at::AbstractAtom) = !iszero(sum(magnetization(at)))

"Takes a Vector of atoms and returns a Vector with the atoms having unique symbols."
function Base.unique(atoms::Vector{A}) where {A<:AbstractAtom}
    uni = A[]
    for at in atoms
        if findfirst(x -> isequal_species(x, at), uni) === nothing
            push!(uni, at)
        end
    end
    return uni
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

isdefault(x::DFTU{T}) where {T} = x == DFTU{T}()

isdefault(x::Any) = isempty(x)
isdefault(x::AbstractVector) = isempty(x) || all(iszero, x)

Base.string(::Type{Elk}, dftu::DFTU) = "$(dftu.l) $(dftu.U) $(dftu.J0)"

Base.isempty(p::Pseudo) = isempty(p.name) && isempty(p.dir)

==(p1::Pseudo, p2::Pseudo) = p1.name == p2.name && p1.dir == p2.dir

path(p::Pseudo) = joinpath(p.dir, p.name)
#Easiest way to implement a new abstractatom is to provide a way to access
#the struct holding `name`, `position_cart`, `element`, `pseudo`, `projection` fields
atom(at::Atom) = at

for interface_function in fieldnames(Atom)
    @eval $interface_function(at::AbstractAtom) = atom(at).$interface_function
    @eval export $interface_function
end

"""
    atoms(job::DFJob)

Returns the atoms inside the structure of the job.
"""
atoms(job::DFJob) = atoms(job.structure)
atoms(job::DFJob, name::Symbol) = job.structure[name]
atoms(job::DFJob, el::Element) = job.structure[el]
atoms(f::Function, job::DFJob) = atoms(f, job.structure)

"""
    atoms([f::Function,], structure::AbstractStructure)
    
Returns `structure.atoms` or `filter(f, structure.atoms)` if `f` is specified.
"""
atoms(str::AbstractStructure) = structure(str).atoms
atoms(f::Function, str::AbstractStructure) = filter(f, atoms(str))
atoms(str::AbstractStructure, el::Element) = str[el]

"""job.structure.atoms = atoms"""
set_atoms!(job::DFJob, atoms::Vector{<:AbstractAtom}) = job.structure.atoms = atoms

function Base.range(at::AbstractAtom)
    projs = projections(at)
    @assert length(projs) != 0 "At $(name(at)) has no defined projections. Please use `setprojections!` first."
    return projs[1].start:projs[end].last
end

Base.range(v::Vector{AbstractAtom}) = vcat(range.(v)...)

set_name!(at::AbstractAtom, name::Symbol) = atom(at).name = name

length_unit(at::Atom{T,LT}) where {T,LT} = LT

function bondlength(at1::AbstractAtom{T}, at2::AbstractAtom{T},
                    R = T(0.0)) where {T<:AbstractFloat}
    return norm(position_cart(at1) - position_cart(at2) - R)
end

function ==(at1::AbstractAtom{T,LT}, at2::AbstractAtom{T,LT}) where {T,LT}
    return name(at1) == name(at2) &&
           norm(position_cart(at1) - position_cart(at2)) < LT(1e-6)
end

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

function position_string(::Type{QE}, at::AbstractAtom; relative = true)
    pos = relative ? round.(position_cryst(at), digits = 5) :
          round.(ustrip.(uconvert.(Ang, position_cart(at))), digits = 5)
    return "$(name(at))  $(pos[1]) $(pos[2]) $(pos[3])\n"
end

"""
    set_magnetization!(at::AbstractAtom, mag; print=true)

Sets the magnetization of the [`Atom`](@ref Atom).

    set_magnetization!(str::Structure, atsym_mag::Pair{Union{Symbol,Element},<:AbstractVector}...)
    set_magnetization!(job::DFJob, atsym_mag::Pair{Union{Symbol,Element},<:AbstractVector}...)

Each of the names in `atsym_mag` will be matched with the [`Atoms`](@ref Atom) and their magnetization
will be set to the specified ones.

Example:
```
set_magnetization!(job, :Ni1 => [0.0, 0.0, -1.0], :Ni2 => [0.0, 0.0, 1.0])
set_magnetization!(job, element(:Ni) => [0.0, 0.0, -1.0])
```
The former will set all the moments of [`Atoms`](@ref Atom) with name `Ni1` to `[0.0, 0.0, -1.0]` and `Ni2` to `[0.0, 0.0, 1.0]`.
The latter will result in all aligned ferromagnetic moments for the same calculation.
Since the moments are aligned with the z-direction this will signal that colinear calculations should be ran, which will
be used when generating the input files.
"""
function set_magnetization!(at::AbstractAtom, mag; print = true)
    at.magnetization = convert(Vec3, mag)
    return print && @info "Magnetization of at $(name(at)) was set to $(magnetization(at))"
end

function set_magnetization!(str::Structure, atsym_mag::Pair{Symbol,<:AbstractVector}...; kwargs...)
    for (atsym, mag) in atsym_mag
        for at in str[atsym]
            set_magnetization!(at, mag; kwargs...)
        end
    end
end

set_magnetization!(job::DFJob, args...) = set_magnetization!(job.structure, args...)

projections(str::AbstractStructure) = projections.(atoms(str))
hasprojections(str::AbstractStructure) = !all(isempty, projections(str))
function hasprojections_assert(str::AbstractStructure)
    @assert hasprojections(str) "No projections found in structure $(str.name).
    Please set the projections for the atoms inside the structure first using `setprojections!`."
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
function set_position!(at::AbstractAtom, pos::AbstractVector{T},
                       unit_cell::Mat3) where {T<:Real}
    atom(at).position_cryst = Point3{T}(pos...)
    atom(at).position_cart  = unit_cell * atom(at).position_cryst
    return at
end

function set_position!(at::AbstractAtom, pos::AbstractVector{T},
                       unit_cell::Mat3) where {T<:Length}
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
    mid = (p1 + p2) / 2
    orig_dist = norm(p1 - p2)
    direction = (p1 - mid) / norm(p1 - mid)
    new_dist = orig_dist * scale
    new_p1 = mid + new_dist * direction / 2
    new_p2 = mid - new_dist * direction / 2
    set_position!(at1, new_p1, cell)
    return set_position!(at2, new_p2, cell)
end

for hub_param in (:U, :J0, :α, :β)
    f = Symbol("set_Hubbard_$(hub_param)!")
    str = "$hub_param"
    @eval begin
        """
            $($(f))(at::AbstractAtom, v::AbstractFloat; print=true)
            $($(f))(str::Structure, ats_$($(str))::Pair{<:Union{Symbol,Element}, <:AbstractFloat}...; print=true)
            $($(f))(job::DFJob, ats_$($(str))::Pair{<:Union{Symbol,Element}, <:AbstractFloat}...; print=true)

        Set the Hubbard $($(str)) parameter for the specified [Atoms](@ref Atom).
        The latter function allows for conveniently setting the parameter for all
        [Atoms](@ref Atom) with the specified `name`.

        Example:
            `$($(f))(job, :Ir => 2.1, element(:Ni) => 1.0, :O => 0.0)`
        """
        function $f(at::AbstractAtom{T}, v::AbstractFloat; print = true) where {T}
            dftu(at).$(hub_param) = convert(T, v)
            return print && @info "Hubbard $($(str)) of atom $(at.name) set to $v"
        end
        function $f(str::Union{DFJob, Structure}, $(hub_param)::Pair{<:Union{Symbol,Element},<:AbstractFloat}...; print = true)
            for (atsym, val) in $(hub_param)
                $f.(atoms(str, atsym), val; print = print)
            end
        end
        export $f
    end
end

"""
	set_Hubbard_J!(at::AbstractAtom, v::Vector{<:AbstractFloat}; print=true)
    set_Hubbard_J!(job::DFJob, ats_Js::Pair{<:Union{Symbol,Element}, Vector{<:AbstractFloat}}...; print=true)
    
Set the Hubbard J parameter for the specified atom.
The latter function allows for conveniently setting the `Hubbard_J` for all
[Atoms](@ref Atom) with the specified `name`.

Example:
	`set_Hubbard_J(at, [2.1])'
    `set_Hubbard_J(job, :Ir => [2.1], :Ni => [1.0])'
"""
function set_Hubbard_J!(at::AbstractAtom{T}, v::Vector{<:AbstractFloat};
                        print = true) where {T}
    dftu(at).J = convert.(T, v)
    return print && @info "Hubbard J of atom $(at.name) set to $v"
end
function set_Hubbard_J!(str::Union{DFJob, Structure}, ats_Js::Pair{<:Union{Symbol,Element}, Vector{<:AbstractFloat}}...;
                        print = true)
    for (atsym, val) in ats_Js
        set_Hubbard_J!.(atoms(str, atsym), val; print = print)
    end
end

export set_Hubbard_J!
