
abstract type AbstractStructure{T, LT} end

mutable struct Structure{T <: AbstractFloat, AA<:AbstractAtom{T}, LT <: Length{T}} <: AbstractStructure{T, LT}
    name ::AbstractString
    cell ::Mat3{LT}
    atoms::Vector{AA}
    data ::Dict{Symbol, Any}
end
Structure() =
	Structure("NoName", eye(3), Atom[], Dict{Symbol, Any}())

Structure(name, cell::Mat3{LT}, atoms::Vector{Atom{T, LT}}) where {T<:AbstractFloat,LT<:Length{T}} =
	Structure{T, Atom{T, LT}, LT}(name, cell, atoms, Dict{Symbol, Any}())

Structure(str::AbstractStructure{T}, atoms::Vector{AT}) where {T<:AbstractFloat, AT<:AbstractAtom{T}} =
	Structure{T, AT}(name(str), cell(str), atoms, data(str))

Structure(cell::Matrix{LT}, atoms::Vector{Atom{T, LT}}) where {T<:AbstractFloat,LT<:Length{T}} =
	Structure{T, Atom{T, LT}, LT}("NoName", cell, atoms, Dict{Symbol, Any}())

Structure(cif_file::String; name="NoName") =
	cif2structure(cif_file, structure_name = name)

structure(str::Structure) = str
"""
Returns all the atoms inside the structure with the specified symbol
"""
function atoms(str::AbstractStructure, atsym::Symbol)
    out = AbstractAtom[]
    for at in str.atoms
        name(at) == atsym && push!(out, at)
    end
    return out
end
atoms(str::AbstractStructure) = structure(str).atoms
name(str::AbstractStructure) = structure(str).name
data(str::AbstractStructure) = structure(str).data

Base.length(str::AbstractStructure) = length(atoms(str))
cell(str::AbstractStructure) = structure(str).cell

projections(str::AbstractStructure) = projections.(atoms(str))
hasprojections(str::AbstractStructure) = !all(isempty, projections(str))
hasprojections_assert(str::AbstractStructure) =
    @assert hasprojections(str) "No projections found in structure $(str.name).
    Please set the projections for the atoms inside the structure first using `setprojections!`."

reciprocal(cell::AbstractMatrix) = inv(cell)

"""
    setprojections!(str::Structure, projs::Pair...)

Sets the projections of the specified atoms. `projs` has to be of form `:atsym => [:proj]`,
where proj = :s, :p, :d, :f, etc.
"""
function setprojections!(str::Structure, projs::Pair...)
    projdict = Dict(projs)
    for at in unique(str.atoms)
        if !haskey(projdict, name(at))
            projdict[name(at)] = [proj.orb for proj in projections(at)]
        end
    end
    emptyprojections!(str)
    addprojections!(atoms(str), projdict)
end

function emptyprojections!(str::Structure)
    for at in str.atoms
        empty!(projections(at))
    end
end

function nprojections(structure)
    n = 0
    for at in atoms(structure)
        projs = projections(at)
        if !isempty(projs)
            n += sum(orbsize.(projs))
        end
    end
    return n
end

#TODO extend base.merge
"Takes a vector of structures and merges all the attributes of the atoms."
function mergestructures(structures::Vector{<:AbstractStructure})
    nonvoid = filter(x -> x != nothing, structures)
    out = nonvoid[1]
    for structure in nonvoid[2:end]
        for at1 in atoms(out), at2 in atoms(structure)
            if at1 == at2
                for fname in fieldnames(typeof(at1))
                    if fname in [:name, :element, :position_cart, :position_cryst]
                        continue
                    end
                    field = getfield(at2, fname)
                    if field == nothing || isdefault(field)
                        continue
                    else
                        setfield!(at1, fname, field)
                    end
                end
            end
        end
    end
    return out
end

"Uses cif2cell to Meta.parse a cif file, then returns the parsed structure."
function cif2structure(cif_file::String; structure_name="NoName")
    tmpdir = dirname(cif_file)
    tmpfile = joinpath(tmpdir, "tmp.in")
    @assert splitext(cif_file)[2] == ".cif" error("Please specify a valid cif input file")
    run(`$pythonpath $cif2cellpath $cif_file -p quantum-espresso -o $tmpfile`)

    bla, structure = qe_read_input(tmpfile, structure_name = structure_name)
    rm(tmpfile)
    return structure
end

function setpseudos!(structure::AbstractStructure, atoms::Vector{<:AbstractAtom}, set::Symbol, specifier::String=""; kwargs...)
	dir = getdefault_pseudodir(set)
	if dir == nothing
		@warn "No pseudos found for set $set."
		return
	end
    for (i, at) in enumerate(atoms)
        pseudo = getdefault_pseudo(name(at), set, specifier=specifier)
        if pseudo == nothing
            @warn "Pseudo for $(name(at)) at index $i not found in set $set."
        else
            setpseudo!(at, pseudo; kwargs...)
        end
    end
end

setpseudos!(structure::AbstractStructure, atname::Symbol, set::Symbol, specifier::String=""; kwargs...) =
	setpseudos!(structure, atoms(structure, atname), set, specifier; kwargs...)

setpseudos!(structure::AbstractStructure, set::Symbol, specifier::String=""; kwargs...) =
	setpseudos!(structure, atoms(structure), set, specifier; kwargs...)

function setpseudos!(structure::AbstractStructure, at_pseudos::Pair{Symbol, Pseudo}...; kwargs...)
    for (atsym, pseudo) in at_pseudos
        for at in atoms(structure, atsym)
            setpseudo!(at, pseudo; kwargs...)
        end
    end
end

"""
    create_supercell(structure::AbstractStructure, na::Int, nb::Int, nc::Int;make_afm=false)

Takes a structure and creates a supercell from it with the given amount of extra cells (`na, nb, nc`) along the a, b, c direction.
If `make_afm` is set to `true` all the labels and magnetizations of the magnetic atoms will be reversed.
"""
function create_supercell(structure::AbstractStructure, na::Int, nb::Int, nc::Int; make_afm=false)
    orig_ats   = atoms(structure)
    orig_cell  = cell(structure)
    scale_mat  = diagm(0 => 1 .+ [na, nb, nc])

    orig_uats = unique(orig_ats)

    new_cell   = orig_cell * scale_mat
    new_atoms  = eltype(orig_ats)[]
    for ia=0:na, ib=0:nb, ic=0:nc
        transl_vec = orig_cell * [ia, ib, ic]
        factor = isodd(ia + ib + ic) ? -1 : 1
        for at in orig_ats
	        cart_pos = position_cart(at) + transl_vec
	        cryst_pos = inv(new_cell) * cart_pos
            if make_afm && ismagnetic(at)
                new_magnetization = factor*magnetization(at)
                push!(new_atoms, Atom(at, projections=deepcopy(projections(at)), magnetization=new_magnetization, position_cart = cart_pos, position_cryst = Point3(cryst_pos)))
            else
                push!(new_atoms, Atom(at, projections=deepcopy(projections(at)), position_cart=cart_pos, position_cryst=Point3(cryst_pos)))
            end
        end
    end
    return Structure(name(structure), Mat3(new_cell), new_atoms, data(structure))
end

"Rescales the cell of the structure."
function scale_cell!(structure::Structure, v)
	scalemat = [v 0 0; 0 v 0; 0 0 v]
	structure.cell *= scalemat
	for at in atoms(structure)
		at.position_cart = structure.cell * at.position_cryst
	end
	structure.cell
end

function set_magnetization!(str::Structure, atsym_mag::Pair{Symbol, <:AbstractVector}...)
	for (atsym, mag) in atsym_mag
		for at in atoms(str, atsym)
			set_magnetization!(at, mag)
		end
	end
end

"""
    volume(cell::Mat3)
	volume(str::Structure)

Calculates the volume for the unit cell.
"""
volume(cell::Mat3) = det(cell)
   
volume(str::Structure) = cell(str)

const DEFAULT_TOLERANCE = 1e-5

struct SPGStructure
    lattice::Matrix{Cdouble}
    positions::Matrix{Cdouble}
    species_indices::Vector{Cint}
end

function SPGStructure(s::Structure)
    clattice = convert(Matrix{Cdouble}, ustrip.(cell(s)'))
    cpositions = convert(Matrix{Cdouble}, hcat(position_cryst.(atoms(s))...))
    uats         = unique(atoms(s))
    species_indices = Cint[findfirst(x -> isequal_species(x, at), uats) for at in atoms(s)]
    return SPGStructure(clattice, cpositions, species_indices)
end

function symmetry_operators(s::SPGStructure; maxsize=52, tolerance = DEFAULT_TOLERANCE)
    rotations    = Array{Cint}(undef, 3, 3, maxsize)
    translations = Array{Cdouble}(undef, 3, maxsize)

    num_ops = ccall((:spg_get_symmetry, SPGLIB), Cint,
      (Ptr{Cint}, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Cint, Cdouble),
       rotations, translations, maxsize, s.lattice, s.positions, s.species_indices, length(s.species_indices), tolerance)
    return (rotations=[Mat3{Int}(rotations[:, :, i]) for i=1:num_ops], translations=[Vec3(translations[:, i]) for i=1:num_ops])
end
symmetry_operators(s::Structure; kwargs...) = symmetry_operators(SPGStructure(s); kwargs...)

function international(s::SPGStructure; tolerance = DEFAULT_TOLERANCE)
    res = zeros(Cchar, 11)

    num = ccall((:spg_get_international, SPGLIB), Cint,
                   (Ptr{Cchar}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Cint, Cdouble),
                   res, s.lattice, s.positions, s.species_indices, length(s.species_indices), tolerance)
    num == 0 && error("Could not determine the international symbol.")

    return num, join(convert(Vector{Char}, res[1:findfirst(iszero, res) - 1]))
end
international(s::Structure; kwargs...) = international(SPGStructure(s); kwargs...)

function niggli_reduce(c::Mat3{Float64}; tolerance = DEFAULT_TOLERANCE)
    ccell = convert(Matrix{Cdouble}, c)
    numops = ccall((:spg_niggli_reduce, SPGLIB), Cint,
                   (Ptr{Cdouble}, Cdouble),
                   ccell, tolerance)
    return Mat3{Float64}(ccell)
end

function niggli_reduce!(s::SPGStructure; tolerance = DEFAULT_TOLERANCE)

    numops = ccall((:spg_niggli_reduce, SPGLIB), Cint,
                   (Ptr{Cdouble}, Cdouble),
                   s.lattice, tolerance)
    numops == 0 && error("Could not determine the niggli reduced cell.")

    return s
end
niggli_reduce(s::Structure; kwargs...) = niggli_reduce!(SPGStructure(s); kwargs...).*unit(eltype(cell(s)))

function find_primitive!(s::SPGStructure; tolerance = DEFAULT_TOLERANCE)
    numats = ccall((:spg_find_primitive, SPGLIB), Cint,
                   (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Cint, Cdouble),
                   s.lattice, s.positions, s.species_indices, length(s.species_indices), tolerance)
    numats == 0 && error("Could not find the primitive of the supplied structure.")
    return SPGStructure(s.lattice, s.positions[:, 1:numats], s.species_indices[1:numats])
end

function find_primitive(s::Structure; kwargs...)
    uats = unique(atoms(s))
    spg = find_primitive!(SPGStructure(s))
    new_cell = Mat3{Float64}(spg.lattice)*unit(eltype(cell(s)))
    newats = eltype(uats)[]
    for i = 1:length(spg.species_indices)
        tat = deepcopy(uats[spg.species_indices[i]])
        set_position!(tat, Point3(spg.positions[:, i]), new_cell)
        push!(newats, tat)
    end
      
    return Structure(s.name, new_cell,newats, s.data) 
end


"""
    cell_parameters(cell::Mat3)
    cell_parameters(str::Structure)

Parameters (a, b, c, α, β, γ) of the input cell returned in a named tuple.
"""
function cell_parameters(cell::Mat3)
    G = transpose(cell) * cell
    a, b, c = sqrt.(diag(G))
    α = acos(0.5(G[2, 3] + G[3, 2])/(c * b))
    β = acos(0.5(G[3, 1] + G[1, 3])/(a * c))
    γ = acos(0.5(G[1, 2] + G[2, 1])/(a * b))
    a, b, c, α, β, γ
end

cell_parameters(s::Structure) = cell_parameters(ustrip.(cell(s)))

function crystal_kind(s::Structure; tolerance=DEFAULT_TOLERANCE)
    n, sym = international(s)
    f = (i, j) -> i <= n <= j
    cs = [:triclinic => (1, 2), :monoclinic => (3, 15),
          :orthorhombic => (16, 74), :tetragonal => (75, 142),
          :trigonal => (143, 167), :hexagonal => (168, 194),
          :cubic => (195, 230)]
    for (k, v) in cs
        if f(v...)
            return k
        end
    end
    error("Crystal kind could not be determined")
end

function lattice_kind(s::Structure; tolerance=DEFAULT_TOLERANCE)
    n, sym = international(s)
    kind = crystal_kind(s; tolerance=tolerance)
    if n ∈ (146, 148, 155, 160, 161, 166, 167)
        return :rhombohedral
    elseif kind == :trigonal
        return :hexagonal
    else
        return kind
    end
end

function high_symmetry_kpoints(s::Structure; tolerance=DEFAULT_TOLERANCE)
    n, sym       = international(s)
    kind = lattice_kind(s)
    primitive    = find_primitive(s)

    primcell = ustrip.(cell(primitive))

    a, b, c, α, β, γ = cell_parameters(ustrip.(cell(s)))

    symerror = () -> error("Unexpected value for international symbol: $sym")
    if kind == :cubic
        if occursin("P", sym)
            return cubic()
        elseif occursin("F", sym)
            return fcc()
        elseif occursin("I", sym)
            return bcc()
        else
            symerror() 
        end
    elseif kind == :tetragonal
        if occursin("P", sym)
            return tet()
        elseif occursin("I", sym)
            if c < a
                return bctet1(c, a)
            else
                return bctet2(c, a)
            end
        else
            symerror()
        end
    elseif kind == :orthorhombic
        if occursin("P", sym)
            return orc()
        elseif occursin("F", sym)
            if 1 / a^2 > 1 / b^2 + 1 / c^2
                return orcf1(a, b, c)
            elseif 1 / a^2 < 1 / b^2 + 1 / c^2
                return orcf2(a, b, c)
            else
                return orcf3(a, b, c)
            end

        elseif occursin("I", sym)
            return orci(a, b, c)

        elseif occursin("C", sym) || occursin("A", sym)
            return orcc(a, b, c)
        else
            symerror()
        end

    elseif kind == :hexagonal
        return hex()

    elseif kind == :rhombohedral
        if α < π/2
            return rhl1(α)
        else
            return rhl2(α)
        end

    elseif kind == :monoclinic
        prim_rec = inv(primcell)
        ka, kb, kc, kα, kβ, kγ = cell_parameters(prim_rec)
        if occursin("P", sym)
            return mcl(b, c, α)

        elseif occursin("C", sym)
            if kγ > π/2
                return mclc1(a, b, c, α)
            elseif kγ == π/2
                return mclc2(a, b, c, α)
            elseif kγ < π/2
                if b * cos(α) / c + b^2 * sin(α)^2 / a^2 < 1
                    return mclc3(a, b, c, α)
                elseif b * cos(α) / c + b^2 * sin(α)^2 / a^2 == 1
                    return mclc4(a, b, c, α)
                elseif b * cos(α) / c + b^2 * sin(α)^2 / a^2 > 1
                    return mclc5(a, b, c, α)
                end
            end
        else
            symerror()
        end

    elseif kind == :triclinic
        prim_rec = inv(primcell)
        ka, kb, kc, kα, kβ, kγ = cell_parameters(prim_rec)
        if kα > π/2 && kβ > π/2 && kγ > π/2
            return tria()
        elseif kα < π/2 && kβ < π/2 && kγ < π/2
            return trib()
        elseif kα > π/2 && kβ > π/2 && kγ == π/2
            return tria()
        elseif kα < π/2 && kβ < π/2 && kγ == π/2
            return trib()
        end


    else
        error("Unknown lattice type $kind")
    end
end

function cubic()
    kpoints = Dict(:Γ => Vec3([0.0, 0.0, 0.0]),
                   :X => Vec3([0.0, 0.5, 0.0]),
                   :R => Vec3([0.5, 0.5, 0.5]),
                   :M => Vec3([0.5, 0.5, 0.0]))
    path = [[:Γ, :X, :M, :Γ, :R, :X], [:M, :R]]
    return (kpoints=kpoints, path=path)
end

function fcc()
    kpoints = Dict(:Γ => Vec3([0.0, 0.0, 0.0]),
                   :K => Vec3([3.0 / 8.0, 3.0 / 8.0, 3.0 / 4.0]),
                   :L => Vec3([0.5, 0.5, 0.5]),
                   :U => Vec3([5.0 / 8.0, 1.0 / 4.0, 5.0 / 8.0]),
                   :W => Vec3([0.5, 1.0 / 4.0, 3.0 / 4.0]),
                   :X => Vec3([0.5, 0.0, 0.5]))

    path = [[:Γ, :X, :W, :K,
                 :Γ, :L, :U, :W, :L, :K], [:U, :X]]
    return (kpoints=kpoints, path=path)
end

function bcc()
    kpoints = Dict(:Γ => Vec3([0.0, 0.0, 0.0]),
                   :H => Vec3([0.5, -0.5, 0.5]),
                   :P => Vec3([0.25, 0.25, 0.25]),
                   :N => Vec3([0.0, 0.0, 0.5]))
    path = [[:Γ, :H, :N, :Γ, :P, :H], [:P, :N]]
    return (kpoints=kpoints, path=path)
end

function tet()
    kpoints = Dict(:Γ => Vec3([0.0, 0.0, 0.0]),
                   :A => Vec3([0.5, 0.5, 0.5]),
                   :M => Vec3([0.5, 0.5, 0.0]),
                   :R => Vec3([0.0, 0.5, 0.5]),
                   :X => Vec3([0.0, 0.5, 0.0]),
                   :Z => Vec3([0.0, 0.0, 0.5]))
        path = [[:Γ, :X, :M, :Γ, :Z, :R, :A, :Z], [:X, :R],
                [:M, :A]]

    return (kpoints=kpoints, path=path)
end

function bctet1(c, a)
    eta = (1 + c^2 / a^2) / 4.0
    kpoints = Dict(:Γ => Vec3([0.0, 0.0, 0.0]),
               :M => Vec3([-0.5, 0.5, 0.5]),
               :N => Vec3([0.0, 0.5, 0.0]),
               :P => Vec3([0.25, 0.25, 0.25]),
               :X => Vec3([0.0, 0.0, 0.5]),
               :Z => Vec3([eta, eta, -eta]),
               :Z_1 => Vec3([-eta, 1 - eta, eta]))
    path = [[:Γ, :X, :M, :Γ, :Z, :P, :N, :Z_1, :M],
            [:X, :P]]
    return (kpoints=kpoints, path=path)
end


function bctet2(c, a)
    eta = (1 + a^2 / c^2) / 4.0
    zeta = a^2 / (2 * c^2)
    kpoints = Dict(:Γ => Vec3([0.0, 0.0, 0.0]),
               :N => Vec3([0.0, 0.5, 0.0]),
               :P => Vec3([0.25, 0.25, 0.25]),
               :Sigma => Vec3([-eta, eta, eta]),
               :Sigma_1 => Vec3([eta, 1 - eta, -eta]),
               :X => Vec3([0.0, 0.0, 0.5]),
               :Y => Vec3([-zeta, zeta, 0.5]),
               :Y_1 => Vec3([0.5, 0.5, -zeta]),
               :Z => Vec3([0.5, 0.5, -0.5]))
    path = [[:Γ, :X, :Y, :Sigma, :Γ, :Z,
             :Sigma_1, :N, :P, :Y_1, :Z], [:X, :P]]
    return (kpoints=kpoints, path=path)
end


function orc()
    kpoints = Dict(:Γ => Vec3([0.0, 0.0, 0.0]),
               :R => Vec3([0.5, 0.5, 0.5]),
               :S => Vec3([0.5, 0.5, 0.0]),
               :T => Vec3([0.0, 0.5, 0.5]),
               :U => Vec3([0.5, 0.0, 0.5]),
               :X => Vec3([0.5, 0.0, 0.0]),
               :Y => Vec3([0.0, 0.5, 0.0]),
               :Z => Vec3([0.0, 0.0, 0.5]))
    path = [[:Γ, :X, :S, :Y, :Γ,
             :Z, :U, :R, :T, :Z], [:Y, :T], [:U, :X], [:S, :R]]
    return (kpoints=kpoints, path=path)
end


function orcf1(a, b, c)
    zeta = (1 + a^2 / b^2 - a^2 / c^2) / 4
    eta = (1 + a^2 / b^2 + a^2 / c^2) / 4

    kpoints = Dict(:Γ => Vec3([0.0, 0.0, 0.0]),
               :A => Vec3([0.5, 0.5 + zeta, zeta]),
               :A_1 => Vec3([0.5, 0.5 - zeta, 1 - zeta]),
               :L => Vec3([0.5, 0.5, 0.5]),
               :T => Vec3([1, 0.5, 0.5]),
               :X => Vec3([0.0, eta, eta]),
               :X_1 => Vec3([1, 1 - eta, 1 - eta]),
               :Y => Vec3([0.5, 0.0, 0.5]),
               :Z => Vec3([0.5, 0.5, 0.0]))
    path = [[:Γ, :Y, :T, :Z, :Γ, :X, :A_1, :Y],
            [:T, :X_1], [:X, :A, :Z], [:L, :Γ]]
    return (kpoints=kpoints, path=path)
end


function orcf2(a, b, c)
    phi = (1 + c^2 / b^2 - c^2 / a^2) / 4
    eta = (1 + a^2 / b^2 - a^2 / c^2) / 4
    delta = (1 + b^2 / a^2 - b^2 / c^2) / 4
    kpoints = Dict(:Γ => Vec3([0.0, 0.0, 0.0]),
               :C => Vec3([0.5, 0.5 - eta, 1 - eta]),
               :C_1 => Vec3([0.5, 0.5 + eta, eta]),
               :D => Vec3([0.5 - delta, 0.5, 1 - delta]),
               :D_1 => Vec3([0.5 + delta, 0.5, delta]),
               :L => Vec3([0.5, 0.5, 0.5]),
               :H => Vec3([1 - phi, 0.5 - phi, 0.5]),
               :H_1 => Vec3([phi, 0.5 + phi, 0.5]),
               :X => Vec3([0.0, 0.5, 0.5]),
               :Y => Vec3([0.5, 0.0, 0.5]),
               :Z => Vec3([0.5, 0.5, 0.0]))
    path = [[:Γ, :Y, :C, :D, :X, :Γ,
             :Z, :D_1, :H, :C], [:C_1, :Z], [:X, :H_1], [:H, :Y],
            [:L, :Γ]]
    return (kpoints=kpoints, path=path)
end


function orcf3(a, b, c)
    zeta = (1 + a^2 / b^2 - a^2 / c^2) / 4
    eta = (1 + a^2 / b^2 + a^2 / c^2) / 4
    kpoints = Dict(:Γ => Vec3([0.0, 0.0, 0.0]),
               :A => Vec3([0.5, 0.5 + zeta, zeta]),
               :A_1 => Vec3([0.5, 0.5 - zeta, 1 - zeta]),
               :L => Vec3([0.5, 0.5, 0.5]),
               :T => Vec3([1, 0.5, 0.5]),
               :X => Vec3([0.0, eta, eta]),
               :X_1 => Vec3([1, 1 - eta, 1 - eta]),
               :Y => Vec3([0.5, 0.0, 0.5]),
               :Z => Vec3([0.5, 0.5, 0.0]))
    path = [[:Γ, :Y, :T, :Z, :Γ, :X, :A_1, :Y],
            [:X, :A, :Z], [:L, :Γ]]
    return (kpoints=kpoints, path=path)
end


function orci(a, b, c)
    zeta = (1 + a^2 / c^2) / 4
    eta = (1 + b^2 / c^2) / 4
    delta = (b^2 - a^2) / (4 * c^2)
    mu = (a^2 + b^2) / (4 * c^2)
    kpoints = Dict(:Γ => Vec3([0.0, 0.0, 0.0]),
               :L => Vec3([-mu, mu, 0.5 - delta]),
               :L_1 => Vec3([mu, -mu, 0.5 + delta]),
               :L_2 => Vec3([0.5 - delta, 0.5 + delta, -mu]),
               :R => Vec3([0.0, 0.5, 0.0]),
               :S => Vec3([0.5, 0.0, 0.0]),
               :T => Vec3([0.0, 0.0, 0.5]),
               :W => Vec3([0.25, 0.25, 0.25]),
               :X => Vec3([-zeta, zeta, zeta]),
               :X_1 => Vec3([zeta, 1 - zeta, -zeta]),
               :Y => Vec3([eta, -eta, eta]),
               :Y_1 => Vec3([1 - eta, eta, -eta]),
               :Z => Vec3([0.5, 0.5, -0.5]))
    path = [[:Γ, :X, :L, :T, :W, :R, :X_1, :Z,
             :Γ, :Y, :S, :W], [:L_1, :Y], [:Y_1, :Z]]
    return (kpoints=kpoints, path=path)
end


function orcc(a, b, c)
    zeta = (1 + a^2 / b^2) / 4
    kpoints = Dict(:Γ => Vec3([0.0, 0.0, 0.0]),
               :A => Vec3([zeta, zeta, 0.5]),
               :A_1 => Vec3([-zeta, 1 - zeta, 0.5]),
               :R => Vec3([0.0, 0.5, 0.5]),
               :S => Vec3([0.0, 0.5, 0.0]),
               :T => Vec3([-0.5, 0.5, 0.5]),
               :X => Vec3([zeta, zeta, 0.0]),
               :X_1 => Vec3([-zeta, 1 - zeta, 0.0]),
               :Y => Vec3([-0.5, 0.5, 0]),
               :Z => Vec3([0.0, 0.0, 0.5]))
    path = [[:Γ, :X, :S, :R, :A, :Z,
             :Γ, :Y, :X_1, :A_1, :T, :Y], [:Z, :T]]
    return (kpoints=kpoints, path=path)
end


function hex()
    kpoints = Dict(:Γ => Vec3([0.0, 0.0, 0.0]),
               :A => Vec3([0.0, 0.0, 0.5]),
               :H => Vec3([1.0 / 3.0, 1.0 / 3.0, 0.5]),
               :K => Vec3([1.0 / 3.0, 1.0 / 3.0, 0.0]),
               :L => Vec3([0.5, 0.0, 0.5]),
               :M => Vec3([0.5, 0.0, 0.0]))
    path = [[:Γ, :M, :K, :Γ, :A, :L, :H, :A], [:L, :M],
            [:K, :H]]
    return (kpoints=kpoints, path=path)
end


function rhl1(alpha)
    eta = (1 + 4 * cos(alpha)) / (2 + 4 * cos(alpha))
    nu = 3.0 / 4.0 - eta / 2.0
    kpoints = Dict(:Γ => Vec3([0.0, 0.0, 0.0]),
               :B => Vec3([eta, 0.5, 1.0 - eta]),
               :B_1 => Vec3([1.0 / 2.0, 1.0 - eta, eta - 1.0]),
               :F => Vec3([0.5, 0.5, 0.0]),
               :L => Vec3([0.5, 0.0, 0.0]),
               :L_1 => Vec3([0.0, 0.0, -0.5]),
               :P => Vec3([eta, nu, nu]),
               :P_1 => Vec3([1.0 - nu, 1.0 - nu, 1.0 - eta]),
               :P_2 => Vec3([nu, nu, eta - 1.0]),
               :Q => Vec3([1.0 - nu, nu, 0.0]),
               :X => Vec3([nu, 0.0, -nu]),
               :Z => Vec3([0.5, 0.5, 0.5]))
    path = [[:Γ, :L, :B_1], [:B, :Z, :Γ, :X],
            [:Q, :F, :P_1, :Z], [:L, :P]]
    return (kpoints=kpoints, path=path)
end


function rhl2(alpha)
    eta = 1 / (2 * tan(alpha / 2.0)^2)
    nu = 3.0 / 4.0 - eta / 2.0
    kpoints = Dict(:Γ => Vec3([0.0, 0.0, 0.0]),
               :F => Vec3([0.5, -0.5, 0.0]),
               :L => Vec3([0.5, 0.0, 0.0]),
               :P => Vec3([1 - nu, -nu, 1 - nu]),
               :P_1 => Vec3([nu, nu - 1.0, nu - 1.0]),
               :Q => Vec3([eta, eta, eta]),
               :Q_1 => Vec3([1.0 - eta, -eta, -eta]),
               :Z => Vec3([0.5, -0.5, 0.5]))
    path = [[:Γ, :P, :Z, :Q, :Γ,
             :F, :P_1, :Q_1, :L, :Z]]
    return (kpoints=kpoints, path=path)
end


function mcl(b, c, beta)
    eta = (1 - b * cos(beta) / c) / (2 * sin(beta)^2)
    nu = 0.5 - eta * c * cos(beta) / b
    kpoints = Dict(:Γ => Vec3([0.0, 0.0, 0.0]),
               :A => Vec3([0.5, 0.5, 0.0]),
               :C => Vec3([0.0, 0.5, 0.5]),
               :D => Vec3([0.5, 0.0, 0.5]),
               :D_1 => Vec3([0.5, 0.5, -0.5]),
               :E => Vec3([0.5, 0.5, 0.5]),
               :H => Vec3([0.0, eta, 1.0 - nu]),
               :H_1 => Vec3([0.0, 1.0 - eta, nu]),
               :H_2 => Vec3([0.0, eta, -nu]),
               :M => Vec3([0.5, eta, 1.0 - nu]),
               :M_1 => Vec3([0.5, 1 - eta, nu]),
               :M_2 => Vec3([0.5, 1 - eta, nu]),
               :X => Vec3([0.0, 0.5, 0.0]),
               :Y => Vec3([0.0, 0.0, 0.5]),
               :Y_1 => Vec3([0.0, 0.0, -0.5]),
               :Z => Vec3([0.5, 0.0, 0.0]))
    path = [[:Γ, :Y, :H, :C, :E, :M_1, :A, :X, :H_1],
            [:M, :D, :Z], [:Y, :D]]
    return (kpoints=kpoints, path=path)
end


function mclc1(a, b, c, alpha)
    zeta = (2 - b * cos(alpha) / c) / (4 * sin(alpha)^2)
    eta = 0.5 + 2 * zeta * c * cos(alpha) / b
    psi = 0.75 - a^2 / (4 * b^2 * sin(alpha)^2)
    phi = psi + (0.75 - psi) * b * cos(alpha) / c
    kpoints = Dict(:Γ => Vec3([0.0, 0.0, 0.0]),
               :N => Vec3([0.5, 0.0, 0.0]),
               :N_1 => Vec3([0.0, -0.5, 0.0]),
               :F => Vec3([1 - zeta, 1 - zeta, 1 - eta]),
               :F_1 => Vec3([zeta, zeta, eta]),
               :F_2 => Vec3([-zeta, -zeta, 1 - eta]),
               # :F_3 => Vec3([1 - zeta, -zeta, 1 - eta]),
               :I => Vec3([phi, 1 - phi, 0.5]),
               :I_1 => Vec3([1 - phi, phi - 1, 0.5]),
               :L => Vec3([0.5, 0.5, 0.5]),
               :M => Vec3([0.5, 0.0, 0.5]),
               :X => Vec3([1 - psi, psi - 1, 0.0]),
               :X_1 => Vec3([psi, 1 - psi, 0.0]),
               :X_2 => Vec3([psi - 1, -psi, 0.0]),
               :Y => Vec3([0.5, 0.5, 0.0]),
               :Y_1 => Vec3([-0.5, -0.5, 0.0]),
               :Z => Vec3([0.0, 0.0, 0.5]))
    path = [[:Γ, :Y, :F, :L, :I], [:I_1, :Z, :F_1],
            [:Y, :X_1], [:X, :Γ, :N], [:M, :Γ]]
    return (kpoints=kpoints, path=path)
end


function mclc2(a, b, c, alpha)
    zeta = (2 - b * cos(alpha) / c) / (4 * sin(alpha)^2)
    eta = 0.5 + 2 * zeta * c * cos(alpha) / b
    psi = 0.75 - a^2 / (4 * b^2 * sin(alpha)^2)
    phi = psi + (0.75 - psi) * b * cos(alpha) / c
    kpoints = Dict(:Γ => Vec3([0.0, 0.0, 0.0]),
               :N => Vec3([0.5, 0.0, 0.0]),
               :N_1 => Vec3([0.0, -0.5, 0.0]),
               :F => Vec3([1 - zeta, 1 - zeta, 1 - eta]),
               :F_1 => Vec3([zeta, zeta, eta]),
               :F_2 => Vec3([-zeta, -zeta, 1 - eta]),
               :F_3 => Vec3([1 - zeta, -zeta, 1 - eta]),
               :I => Vec3([phi, 1 - phi, 0.5]),
               :I_1 => Vec3([1 - phi, phi - 1, 0.5]),
               :L => Vec3([0.5, 0.5, 0.5]),
               :M => Vec3([0.5, 0.0, 0.5]),
               :X => Vec3([1 - psi, psi - 1, 0.0]),
               :X_1 => Vec3([psi, 1 - psi, 0.0]),
               :X_2 => Vec3([psi - 1, -psi, 0.0]),
               :Y => Vec3([0.5, 0.5, 0.0]),
               :Y_1 => Vec3([-0.5, -0.5, 0.0]),
               :Z => Vec3([0.0, 0.0, 0.5]))
    path = [[:Γ, :Y, :F, :L, :I], [:I_1, :Z, :F_1],
            [:N, :Γ, :M]]
    return (kpoints=kpoints, path=path)
end


function mclc3(a, b, c, alpha)
    mu = (1 + b^2 / a^2) / 4.0
    delta = b * c * cos(alpha) / (2 * a^2)
    zeta = mu - 0.25 + (1 - b * cos(alpha) / c) / (4 * sin(alpha)^2)
    eta = 0.5 + 2 * zeta * c * cos(alpha) / b
    phi = 1 + zeta - 2 * mu
    psi = eta - 2 * delta
    kpoints = Dict(:Γ => Vec3([0.0, 0.0, 0.0]),
               :F => Vec3([1 - phi, 1 - phi, 1 - psi]),
               :F_1 => Vec3([phi, phi - 1, psi]),
               :F_2 => Vec3([1 - phi, -phi, 1 - psi]),
               :H => Vec3([zeta, zeta, eta]),
               :H_1 => Vec3([1 - zeta, -zeta, 1 - eta]),
               :H_2 => Vec3([-zeta, -zeta, 1 - eta]),
               :I => Vec3([0.5, -0.5, 0.5]),
               :M => Vec3([0.5, 0.0, 0.5]),
               :N => Vec3([0.5, 0.0, 0.0]),
               :N_1 => Vec3([0.0, -0.5, 0.0]),
               :X => Vec3([0.5, -0.5, 0.0]),
               :Y => Vec3([mu, mu, delta]),
               :Y_1 => Vec3([1 - mu, -mu, -delta]),
               :Y_2 => Vec3([-mu, -mu, -delta]),
               :Y_3 => Vec3([mu, mu - 1, delta]),
               :Z => Vec3([0.0, 0.0, 0.5]))
    path = [[:Γ, :Y, :F, :H, :Z, :I, :F_1],
            [:H_1, :Y_1, :X, :Γ, :N], [:M, :Γ]]
    return (kpoints=kpoints, path=path)
end


function mclc4(a, b, c, alpha)
    mu = (1 + b^2 / a^2) / 4.0
    delta = b * c * cos(alpha) / (2 * a^2)
    zeta = mu - 0.25 + (1 - b * cos(alpha) / c) / (4 * sin(alpha)^2)
    eta = 0.5 + 2 * zeta * c * cos(alpha) / b
    phi = 1 + zeta - 2 * mu
    psi = eta - 2 * delta
    kpoints = Dict(:Γ => Vec3([0.0, 0.0, 0.0]),
               :F => Vec3([1 - phi, 1 - phi, 1 - psi]),
               :F_1 => Vec3([phi, phi - 1, psi]),
               :F_2 => Vec3([1 - phi, -phi, 1 - psi]),
               :H => Vec3([zeta, zeta, eta]),
               :H_1 => Vec3([1 - zeta, -zeta, 1 - eta]),
               :H_2 => Vec3([-zeta, -zeta, 1 - eta]),
               :I => Vec3([0.5, -0.5, 0.5]),
               :M => Vec3([0.5, 0.0, 0.5]),
               :N => Vec3([0.5, 0.0, 0.0]),
               :N_1 => Vec3([0.0, -0.5, 0.0]),
               :X => Vec3([0.5, -0.5, 0.0]),
               :Y => Vec3([mu, mu, delta]),
               :Y_1 => Vec3([1 - mu, -mu, -delta]),
               :Y_2 => Vec3([-mu, -mu, -delta]),
               :Y_3 => Vec3([mu, mu - 1, delta]),
               :Z => Vec3([0.0, 0.0, 0.5]))
    path = [[:Γ, :Y, :F, :H, :Z, :I],
            [:H_1, :Y_1, :X, :Γ, :N], [:M, :Γ]]
    return (kpoints=kpoints, path=path)
end


function mclc5(a, b, c, alpha)
    zeta = (b^2 / a^2 + (1 - b * cos(alpha) / c)
            / sin(alpha)^2) / 4
    eta = 0.5 + 2 * zeta * c * cos(alpha) / b
    mu = eta / 2 + b^2 / (4 * a^2) - b * c * cos(alpha) / (2 * a^2)
    nu = 2 * mu - zeta
    rho = 1 - zeta * a^2 / b^2
    omega = (4 * nu - 1 - b^2 * sin(alpha)^2 / a^2) * c / (2 * b * cos(alpha))
    delta = zeta * c * cos(alpha) / b + omega / 2 - 0.25
    kpoints = Dict(:Γ => Vec3([0.0, 0.0, 0.0]),
               :F => Vec3([nu, nu, omega]),
               :F_1 => Vec3([1 - nu, 1 - nu, 1 - omega]),
               :F_2 => Vec3([nu, nu - 1, omega]),
               :H => Vec3([zeta, zeta, eta]),
               :H_1 => Vec3([1 - zeta, -zeta, 1 - eta]),
               :H_2 => Vec3([-zeta, -zeta, 1 - eta]),
               :I => Vec3([rho, 1 - rho, 0.5]),
               :I_1 => Vec3([1 - rho, rho - 1, 0.5]),
               :L => Vec3([0.5, 0.5, 0.5]),
               :M => Vec3([0.5, 0.0, 0.5]),
               :N => Vec3([0.5, 0.0, 0.0]),
               :N_1 => Vec3([0.0, -0.5, 0.0]),
               :X => Vec3([0.5, -0.5, 0.0]),
               :Y => Vec3([mu, mu, delta]),
               :Y_1 => Vec3([1 - mu, -mu, -delta]),
               :Y_2 => Vec3([-mu, -mu, -delta]),
               :Y_3 => Vec3([mu, mu - 1, delta]),
               :Z => Vec3([0.0, 0.0, 0.5]))
    path = [[:Γ, :Y, :F, :L, :I], [:I_1, :Z, :H, :F_1],
            [:H_1, :Y_1, :X, :Γ, :N], [:M, :Γ]]
    return (kpoints=kpoints, path=path)
end


function tria()
    kpoints = Dict(:Γ => Vec3([0.0, 0.0, 0.0]),
               :L => Vec3([0.5, 0.5, 0.0]),
               :M => Vec3([0.0, 0.5, 0.5]),
               :N => Vec3([0.5, 0.0, 0.5]),
               :R => Vec3([0.5, 0.5, 0.5]),
               :X => Vec3([0.5, 0.0, 0.0]),
               :Y => Vec3([0.0, 0.5, 0.0]),
               :Z => Vec3([0.0, 0.0, 0.5]))
    path = [[:X, :Γ, :Y], [:L, :Γ, :Z],
            [:N, :Γ, :M], [:R, :Γ]]
    return (kpoints=kpoints, path=path)
end


function trib()
    kpoints = Dict(:Γ => Vec3([0.0, 0.0, 0.0]),
               :L => Vec3([0.5, -0.5, 0.0]),
               :M => Vec3([0.0, 0.0, 0.5]),
               :N => Vec3([-0.5, -0.5, 0.5]),
               :R => Vec3([0.0, -0.5, 0.5]),
               :X => Vec3([0.0, -0.5, 0.0]),
               :Y => Vec3([0.5, 0.0, 0.0]),
               :Z => Vec3([-0.5, 0.0, 0.5]))
    path = [[:X, :Γ, :Y], [:L, :Γ, :Z],
            [:N, :Γ, :M], [:R, :Γ]]
    return (kpoints=kpoints, path=path)
end
