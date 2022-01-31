"""
    Structure(cell::Mat3, atoms::Vector{Atom})

The structure on which the [`Calculations`](@ref Calculation) will be performed.

    Structure(cif_file::String)

Creates a [`Structure`](@ref) from the supplied cif file.
"""
mutable struct Structure
    cell  :: Mat3{typeof(1.0Ang)}
    atoms :: Vector{Atom}
end
Structure() = Structure(Mat3(fill(1.0Ang, 3, 3)), Atom[])

function Structure(cif_file::String)
    str = cif2structure(cif_file)
    @info "Structure extracted from $cif_file\n\tcell parameters: \n\t a = $((str.cell[:,1]...,))\n\t b = $((str.cell[:,2]...,))\n\t c = $((str.cell[:,3]...,))\n\tnat = $(length(str.atoms))\n\telements = $(unique(getfield.(str.atoms, :name)))"
    return str
end

StructTypes.StructType(::Type{Structure}) = StructTypes.Struct()

function Base.:(==)(str1::Structure, str2::Structure)
    return str1.cell == str2.cell && str1.atoms == str2.atoms
end

"Uses cif2cell to Meta.parse a cif file, then returns the parsed structure."
function cif2structure(cif_file::String)
    tmpdir = dirname(cif_file)
    tmpfile = joinpath(tmpdir, "tmp.in")
    @assert splitext(cif_file)[2] == ".cif" error("Please specify a valid cif calculation file")
    if Sys.which("cif2cell") === nothing 
        run(`$(DFControl.PYTHONPATH) $(DFControl.CIF2CELLPATH) $cif_file --no-reduce -p quantum-espresso -o $tmpfile`)
    else
        run(`cif2cell $cif_file --no-reduce -p quantum-espresso -o $tmpfile`)
    end
    bla, structure = DFC.FileIO.qe_read_calculation(tmpfile)
    rm(tmpfile)
    return structure
end

# TODO extend base.merge
"Takes a vector of structures and merges all the attributes of the atoms."
function mergestructures(structures::Vector{Structure})
    nonvoid = filter(x -> x != nothing, structures)
    out = nonvoid[1]
    default_at = Atom(; name = :Rh, position_cart = zero(Point3{typeof(1.0Ang)}),
                      position_cryst = zero(Point3{Float64}))
    for structure in nonvoid[2:end]
        for at1 in out.atoms, at2 in structure.atoms
            if at1 == at2
                for fname in fieldnames(typeof(at1))
                    if fname in [:name, :element, :position_cart, :position_cryst]
                        continue
                    end
                    field = getfield(at2, fname)
                    if field != getfield(default_at, fname)
                        setfield!(at1, fname, field)
                    end
                end
            end
        end
    end
    return out
end

"""
    update_geometry!(str1::Structure, str2::Structure)

Updates the spatial parameters of the `atoms` and `cell` of the first structure to those found in the second.
"""
function update_geometry!(str1::Structure, str2::Structure)
    str1.cell = copy(str2.cell)
    ats2 = str2.atoms
    for at1 in str1.atoms
        tats2 = filter(y -> y.name == at1.name, ats2)
        id = findmin(map(x -> norm(x.position_cart - at1.position_cart), tats2))[2]
        if id === nothing
            @error "No atom of the species $(at1.name) found in the second structure"
        end
        set_position!(at1, tats2[id].position_cryst, str1.cell)
    end
    return str1
end
Base.length(str::Structure) = length(str.atoms)
Base.iterate(str::Structure, args...) = iterate(str.atoms, args...)

### CELL ###
"""
    set_cell!(structure::Structure, c::Mat3)

Sets the `cell` of the `Structure` to `c`.
`c` will first be converted to the right units to match those
of the current `cell`.
All the absolute coordinates of the atoms will be recalculated
based on the new `cell`. 
"""
function set_cell!(structure::Structure, c::Mat3)
    @warn "Converting cell to the units of the current cell of the structure."
    cell_ = uconvert.(Ang, c)
    for a in structure.atoms
        set_position!(a, a.position_cryst, cell_)
    end
    return structure.cell = cell_
end

"""
    scale_cell!(structure::Structure, scalemat::Matrix)
    
Rescales the cell of the `structure`.
"""
function scale_cell!(structure::Structure, scalemat::Matrix)
    new_cell = Mat3(structure.cell * scalemat)
    return set_cell!(structure, new_cell)
end

"""
    volume(cell::Mat3)
    volume(str::Structure)

Calculates the volume for the unit cell.
"""
volume(cell::Mat3) = det(cell)
volume(str::Structure) = volume(str.cell)
reciprocal(cell::AbstractMatrix) = inv(cell)

"""
    a(str::Structure)

First lattice vector.
"""
a(str::Structure) = str.cell[:, 1]
"""
    b(str::Structure)

Second lattice vector.
"""
b(str::Structure) = str.cell[:, 2]
"""
    c(str::Structure)

Third lattice vector.
"""
c(str::Structure) = str.cell[:, 3]

"""
    create_supercell(structure::Structure, na::Int, nb::Int, nc::Int; make_afm=false)
    create_supercell(structure::Structure, na::UnitRange, nb::UnitRange, nc::UnitRange; make_afm=false)

Takes a structure and creates a supercell from it with: the given amount of additional cells if (`na::Int, nb::Int, nc::Int`) along the a, b, c direction, or amount of cells specified by the ranges i.e. `-1:1, -1:1, -1:1` would create a 3x3x3 supercell.
If `make_afm` is set to `true` all the labels and magnetizations of the magnetic atoms will be reversed in a checkerboard fashion.
"""
function create_supercell(structure::Structure, na::UnitRange, nb::UnitRange, nc::UnitRange;
                          make_afm = false)
    orig_ats  = structure.atoms
    orig_cell = structure.cell
    scale_mat = diagm(0 => length.([na, nb, nc]))

    orig_uats = unique(orig_ats)

    new_cell  = orig_cell * scale_mat
    new_atoms = eltype(orig_ats)[]
    for ia in na, ib in nb, ic in nc
        transl_vec = orig_cell * [ia, ib, ic]
        factor = isodd(ia + ib + ic) ? -1 : 1
        for at in orig_ats
            cart_pos = at.position_cart + transl_vec
            cryst_pos = inv(new_cell) * cart_pos
            if make_afm && norm(at.magnetization) != 0
                new_magnetization = factor * at.magnetization
                push!(new_atoms,
                      Atom(at; projections = deepcopy(at.projections),
                           magnetization = new_magnetization, position_cart = cart_pos,
                           position_cryst = Point3(cryst_pos)))
            else
                push!(new_atoms,
                      Atom(at; projections = deepcopy(at.projections),
                           position_cart = cart_pos, position_cryst = Point3(cryst_pos)))
            end
        end
    end
    out = Structure(Mat3(new_cell), new_atoms)
    sanitize!(out)
    return out
end

function create_supercell(structure::Structure, na::Int, nb::Int, nc::Int; make_afm = false)
    return create_supercell(structure, 0:na, 0:nb, 0:nc; make_afm = make_afm)
end

### Atoms ###
"""
    getindex(structure::Structure, i::Int)
    getindex(structure::Structure, name::Symbol)
    getindex(structure::Structure, el::Element)

Returns the `i`th atom in `structure`, or all atoms with `name` or are of element `el`.
"""
Base.getindex(str::Structure, args...) = str.atoms[args...]
Base.getindex(str::Structure, el::Symbol) = filter(x -> x.name == el, str.atoms)
Base.getindex(str::Structure, el::Element) = filter(x -> x.element == el, str.atoms)
ismagnetic(str::Structure) = any(x -> norm(x.magnetization) > 0, str.atoms)

## Magnetism ##
function iscolin(str::Structure)
    magats = filter(x -> norm(x.magnetization) > 0, str.atoms)
    return !isempty(magats) &&
           all(x -> isapprox(abs(normalize(x.magnetization) ⋅ [0, 0, 1]), 1.0), magats)
end

function isnoncolin(str::Structure)
    return any(x -> norm(x.magnetization) > 0 &&
                   !isapprox(abs(normalize(x.magnetization) ⋅ [0, 0, 1]), 1.0), str.atoms)
end

function sanitize!(str::Structure)
    magnetic_ats = filter(a -> norm(a.magnetization) != 0 || a.dftu.U != 0.0, str.atoms)
    magnetic_elements = map(x -> x.element, magnetic_ats)
    for e in magnetic_elements
        magnetizations = Vec3[]
        dftus          = DFTU[]
        names          = Symbol[]
        for a in filter(x -> x.element == e, magnetic_ats)
            nameid = findfirst(x -> x == a.name, names)
            if nameid === nothing
                push!(names, a.name)
                push!(magnetizations, a.magnetization)
                push!(dftus, a.dftu)
            else
                mag = magnetizations[nameid]
                dftu_ = dftus[nameid]
                if (a.magnetization != mag || a.dftu != dftu_)
                    id = findfirst(x -> dftus[x] == a.dftu &&
                                       magnetizations[x] == a.magnetization,
                                   1:length(dftus))
                    if id === nothing
                        id = 1
                        tname = Symbol(string(e.symbol) * "$id")
                        while tname ∈ names
                            id += 1
                            tname = Symbol(string(e.symbol) * "$id")
                        end
                        push!(names, tname)
                        push!(magnetizations, a.magnetization)
                        push!(dftus, a.dftu)
                        oldname = a.name
                        a.name = tname
                        @info "Renamed atom from $oldname to $(a.name) in order to distinguish different magnetization species."
                    else
                        oldname = a.name
                        a.name = names[id]
                        @info "Renamed atom from $oldname to $(a.name) in order to distinguish different magnetization species."
                    end
                end
            end
        end
    end
end

"""
    polyhedron(at::Atom, atoms::Vector{Atom}, order::Int)
    polyhedron(at::Atom, str::Structure, order::Int)

Returns a polyhedron around the atom, i.e. the `order` closest atoms.
The returned atoms will be ordered according to their distance to the first one.
In the case of a structure rather than a set of atoms, the search will
be performed over all atoms in the structure.
"""
function polyhedron(at::Atom, atoms::Vector{Atom}, order::Int)
    return sort(atoms; by = x -> norm(x.position_cart - at.position_cart))[1:order]
end
function polyhedron(at::Atom, str::Structure, order::Int)
    return polyhedron(at, create_supercell(str, -1:1, -1:1, -1:1).atoms, order)
end

function set_pseudos!(structure::Structure, pseudos::Dict)
    for at in structure.atoms
        pseudo = get(pseudos, at.element.symbol, nothing)
        if pseudo === nothing
            @warn "Pseudo for $(at.name) not found."
        else
            at.pseudo = pseudo
        end
    end
end

### SPGLIB ###
const DEFAULT_TOLERANCE = 1e-5

struct SPGStructure
    lattice::Matrix{Cdouble}
    positions::Matrix{Cdouble}
    species_indices::Vector{Cint}
end

function SPGStructure(s::Structure)
    clattice = convert(Matrix{Cdouble}, ustrip.(s.cell'))
    cpositions = convert(Matrix{Cdouble}, hcat(map(x -> x.position_cryst, s.atoms)...))
    uats = unique(s.atoms)
    species_indices = Cint[findfirst(x -> isequal_species(x, at), uats) for at in s.atoms]
    return SPGStructure(clattice, cpositions, species_indices)
end

"""
    symmetry_operators(s::Structure; maxsize=52, tolerance=$DEFAULT_TOLERANCE)

Finds and returns all the rotations and translations that are symmetry operators of the structure.
"""
symmetry_operators(s::Structure; kwargs...) = symmetry_operators(SPGStructure(s); kwargs...)
function symmetry_operators(s::SPGStructure; maxsize = 52, tolerance = DEFAULT_TOLERANCE)
    rotations    = Array{Cint}(undef, 3, 3, maxsize)
    translations = Array{Cdouble}(undef, 3, maxsize)

    num_ops = ccall((:spg_get_symmetry, SPGLIB), Cint,
                    (Ptr{Cint}, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint},
                     Cint, Cdouble), rotations, translations, maxsize, s.lattice,
                    s.positions, s.species_indices, length(s.species_indices), tolerance)
    return (rotations = [Mat3{Int}(rotations[:, :, i]) for i in 1:num_ops],
            translations = [Vec3(translations[:, i]) for i in 1:num_ops])
end

"""
    international(s::Structure; tolerance=$DEFAULT_TOLERANCE)

Returns the international symbol of the space group of the structure.
"""
international(s::Structure; kwargs...) = international(SPGStructure(s); kwargs...)
function international(s::SPGStructure; tolerance = DEFAULT_TOLERANCE)
    res = zeros(Cchar, 11)

    num = ccall((:spg_get_international, SPGLIB), Cint,
                (Ptr{Cchar}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Cint, Cdouble), res,
                s.lattice, s.positions, s.species_indices, length(s.species_indices),
                tolerance)
    num == 0 && error("Could not determine the international symbol.")

    return num, join(convert(Vector{Char}, res[1:findfirst(iszero, res)-1]))
end

"""
    niggli_reduce(s::Structure; tolerance=$DEFAULT_TOLERANCE)

Returns the niggli reduced `Structure`.
"""
function niggli_reduce(s::Structure; tolerance = DEFAULT_TOLERANCE)
    reduced_spg_structure = niggli_reduce!(SPGStructure(s); tolerance = DEFAULT_TOLERANCE)
    uats = unique(s.atoms)
    c = Mat3{Float64}(reduced_spg_structure.lattice') .* 1Ang
    return Structure(c, map(1:length(reduced_spg_structure.species_indices)) do i
                         at = deepcopy(uats[reduced_spg_structure.species_indices[i]])
                         pos = reduced_spg_structure.positions[:, i]
                         set_position!(at, pos, c)
                         return at
                     end)
end

function niggli_reduce(c::Mat3{Float64}; tolerance = DEFAULT_TOLERANCE)
    ccell = convert(Matrix{Cdouble}, c)
    numops = ccall((:spg_niggli_reduce, SPGLIB), Cint, (Ptr{Cdouble}, Cdouble), ccell,
                   tolerance)
    return Mat3{Float64}(ccell)
end

function niggli_reduce!(s::SPGStructure; tolerance = DEFAULT_TOLERANCE)
    numops = ccall((:spg_niggli_reduce, SPGLIB), Cint, (Ptr{Cdouble}, Cdouble), s.lattice,
                   tolerance)
    numops == 0 && error("Could not determine the niggli reduced cell.")

    return s
end

function find_primitive!(s::SPGStructure; tolerance = DEFAULT_TOLERANCE)
    numats = ccall((:spg_find_primitive, SPGLIB), Cint,
                   (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Cint, Cdouble), s.lattice,
                   s.positions, s.species_indices, length(s.species_indices), tolerance)
    numats == 0 && error("Could not find the primitive of the supplied structure.")
    return SPGStructure(s.lattice, s.positions[:, 1:numats], s.species_indices[1:numats])
end

function find_primitive(s::Structure; kwargs...)
    uats = unique(s.atoms)
    spg = find_primitive!(SPGStructure(s))
    new_cell = Mat3{Float64}(spg.lattice) * unit(eltype(s.cell))
    newats = eltype(uats)[]
    for i in 1:length(spg.species_indices)
        tat = deepcopy(uats[spg.species_indices[i]])
        set_position!(tat, Point3(spg.positions[:, i]), new_cell)
        push!(newats, tat)
    end

    return Structure(new_cell, newats)
end

"""
    cell_parameters(cell::Mat3)
    cell_parameters(str::Structure)

Parameters `(a, b, c, α, β, γ)`of the calculation cell returned in a `NamedTuple`.
"""
cell_parameters(s::Structure) = cell_parameters(ustrip.(s.cell))
function cell_parameters(cell::Mat3)
    G = transpose(cell) * cell
    a, b, c = sqrt.(diag(G))
    α = acos(0.5(G[2, 3] + G[3, 2]) / (c * b))
    β = acos(0.5(G[3, 1] + G[1, 3]) / (a * c))
    γ = acos(0.5(G[1, 2] + G[2, 1]) / (a * b))
    return a, b, c, α, β, γ
end

function crystal_kind(s::Structure; tolerance = DEFAULT_TOLERANCE)
    n, sym = international(s)
    f = (i, j) -> i <= n <= j
    cs = [:triclinic => (1, 2), :monoclinic => (3, 15), :orthorhombic => (16, 74),
          :tetragonal => (75, 142), :trigonal => (143, 167), :hexagonal => (168, 194),
          :cubic => (195, 230)]
    for (k, v) in cs
        if f(v...)
            return k
        end
    end
    return error("Crystal kind could not be determined")
end

function lattice_kind(s::Structure; tolerance = DEFAULT_TOLERANCE)
    n, sym = international(s)
    kind = crystal_kind(s; tolerance = tolerance)
    if n ∈ (146, 148, 155, 160, 161, 166, 167)
        return :rhombohedral
    elseif kind == :trigonal
        return :hexagonal
    else
        return kind
    end
end

"""
    high_symmetry_kpath(s::Structure, npoints_per_segment::Int; package=QE, kwargs...)

Generates a QE bands calculation compliant high symmetry kpath, to be used with
e.g. `set_kpoints!(bands_calculation, kpoints)`. 
"""
function high_symmetry_kpath(s::Structure, npoints_per_segment::Int; package = QE,
                             kwargs...)
    kpoints, kpath = high_symmetry_kpoints(s)
    out_k = NTuple{4,Float64}[]
    if package === QE
        for path in kpath
            for p in path
                push!(out_k, (kpoints[p]..., float(npoints_per_segment)))
            end
        end
    else
        error("high_symmetry_kpath not implemented for package $package.")
    end
    out_k[end] = (out_k[end][1:3]..., 1.0)
    return out_k
end

"""
    high_symmetry_kpoints(s::Structure; tolerance = 1e-5)

Returns `(kpoints, path)` where `kpoints` are the high-symmetry k-points,
and `path` are the sections of the high symmetry path through the first Brillouin Zone.
"""
function high_symmetry_kpoints(s::Structure; tolerance = DEFAULT_TOLERANCE)
    n, sym    = international(s)
    kind      = lattice_kind(s)
    primitive = find_primitive(s)

    primcell = ustrip.(primitive.cell)

    a, b, c, α, β, γ = cell_parameters(s)

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
        if α < π / 2
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
            if kγ > π / 2
                return mclc1(a, b, c, α)
            elseif kγ == π / 2
                return mclc2(a, b, c, α)
            elseif kγ < π / 2
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
        if kα > π / 2 && kβ > π / 2 && kγ > π / 2
            return tria()
        elseif kα < π / 2 && kβ < π / 2 && kγ < π / 2
            return trib()
        elseif kα > π / 2 && kβ > π / 2 && kγ == π / 2
            return tria()
        elseif kα < π / 2 && kβ < π / 2 && kγ == π / 2
            return trib()
        end

    else
        error("Unknown lattice type $kind")
    end
end

function cubic()
    kpoints = Dict(:Γ => Vec3([0.0, 0.0, 0.0]), :X => Vec3([0.0, 0.5, 0.0]),
                   :R => Vec3([0.5, 0.5, 0.5]), :M => Vec3([0.5, 0.5, 0.0]))
    path = [[:Γ, :X, :M, :Γ, :R, :X], [:M, :R]]
    return (kpoints = kpoints, path = path)
end

function fcc()
    kpoints = Dict(:Γ => Vec3([0.0, 0.0, 0.0]),
                   :K => Vec3([3.0 / 8.0, 3.0 / 8.0, 3.0 / 4.0]),
                   :L => Vec3([0.5, 0.5, 0.5]),
                   :U => Vec3([5.0 / 8.0, 1.0 / 4.0, 5.0 / 8.0]),
                   :W => Vec3([0.5, 1.0 / 4.0, 3.0 / 4.0]), :X => Vec3([0.5, 0.0, 0.5]))

    path = [[:Γ, :X, :W, :K, :Γ, :L, :U, :W, :L, :K], [:U, :X]]
    return (kpoints = kpoints, path = path)
end

function bcc()
    kpoints = Dict(:Γ => Vec3([0.0, 0.0, 0.0]), :H => Vec3([0.5, -0.5, 0.5]),
                   :P => Vec3([0.25, 0.25, 0.25]), :N => Vec3([0.0, 0.0, 0.5]))
    path = [[:Γ, :H, :N, :Γ, :P, :H], [:P, :N]]
    return (kpoints = kpoints, path = path)
end

function tet()
    kpoints = Dict(:Γ => Vec3([0.0, 0.0, 0.0]), :A => Vec3([0.5, 0.5, 0.5]),
                   :M => Vec3([0.5, 0.5, 0.0]), :R => Vec3([0.0, 0.5, 0.5]),
                   :X => Vec3([0.0, 0.5, 0.0]), :Z => Vec3([0.0, 0.0, 0.5]))
    path = [[:Γ, :X, :M, :Γ, :Z, :R, :A, :Z], [:X, :R], [:M, :A]]

    return (kpoints = kpoints, path = path)
end

function bctet1(c, a)
    eta = (1 + c^2 / a^2) / 4.0
    kpoints = Dict(:Γ => Vec3([0.0, 0.0, 0.0]), :M => Vec3([-0.5, 0.5, 0.5]),
                   :N => Vec3([0.0, 0.5, 0.0]), :P => Vec3([0.25, 0.25, 0.25]),
                   :X => Vec3([0.0, 0.0, 0.5]), :Z => Vec3([eta, eta, -eta]),
                   :Z_1 => Vec3([-eta, 1 - eta, eta]))
    path = [[:Γ, :X, :M, :Γ, :Z, :P, :N, :Z_1, :M], [:X, :P]]
    return (kpoints = kpoints, path = path)
end

function bctet2(c, a)
    eta = (1 + a^2 / c^2) / 4.0
    zeta = a^2 / (2 * c^2)
    kpoints = Dict(:Γ => Vec3([0.0, 0.0, 0.0]), :N => Vec3([0.0, 0.5, 0.0]),
                   :P => Vec3([0.25, 0.25, 0.25]), :Sigma => Vec3([-eta, eta, eta]),
                   :Sigma_1 => Vec3([eta, 1 - eta, -eta]), :X => Vec3([0.0, 0.0, 0.5]),
                   :Y => Vec3([-zeta, zeta, 0.5]), :Y_1 => Vec3([0.5, 0.5, -zeta]),
                   :Z => Vec3([0.5, 0.5, -0.5]))
    path = [[:Γ, :X, :Y, :Sigma, :Γ, :Z, :Sigma_1, :N, :P, :Y_1, :Z], [:X, :P]]
    return (kpoints = kpoints, path = path)
end

function orc()
    kpoints = Dict(:Γ => Vec3([0.0, 0.0, 0.0]), :R => Vec3([0.5, 0.5, 0.5]),
                   :S => Vec3([0.5, 0.5, 0.0]), :T => Vec3([0.0, 0.5, 0.5]),
                   :U => Vec3([0.5, 0.0, 0.5]), :X => Vec3([0.5, 0.0, 0.0]),
                   :Y => Vec3([0.0, 0.5, 0.0]), :Z => Vec3([0.0, 0.0, 0.5]))
    path = [[:Γ, :X, :S, :Y, :Γ, :Z, :U, :R, :T, :Z], [:Y, :T], [:U, :X], [:S, :R]]
    return (kpoints = kpoints, path = path)
end

function orcf1(a, b, c)
    zeta = (1 + a^2 / b^2 - a^2 / c^2) / 4
    eta = (1 + a^2 / b^2 + a^2 / c^2) / 4

    kpoints = Dict(:Γ => Vec3([0.0, 0.0, 0.0]), :A => Vec3([0.5, 0.5 + zeta, zeta]),
                   :A_1 => Vec3([0.5, 0.5 - zeta, 1 - zeta]), :L => Vec3([0.5, 0.5, 0.5]),
                   :T => Vec3([1, 0.5, 0.5]), :X => Vec3([0.0, eta, eta]),
                   :X_1 => Vec3([1, 1 - eta, 1 - eta]), :Y => Vec3([0.5, 0.0, 0.5]),
                   :Z => Vec3([0.5, 0.5, 0.0]))
    path = [[:Γ, :Y, :T, :Z, :Γ, :X, :A_1, :Y], [:T, :X_1], [:X, :A, :Z], [:L, :Γ]]
    return (kpoints = kpoints, path = path)
end

function orcf2(a, b, c)
    phi = (1 + c^2 / b^2 - c^2 / a^2) / 4
    eta = (1 + a^2 / b^2 - a^2 / c^2) / 4
    delta = (1 + b^2 / a^2 - b^2 / c^2) / 4
    kpoints = Dict(:Γ => Vec3([0.0, 0.0, 0.0]), :C => Vec3([0.5, 0.5 - eta, 1 - eta]),
                   :C_1 => Vec3([0.5, 0.5 + eta, eta]),
                   :D => Vec3([0.5 - delta, 0.5, 1 - delta]),
                   :D_1 => Vec3([0.5 + delta, 0.5, delta]), :L => Vec3([0.5, 0.5, 0.5]),
                   :H => Vec3([1 - phi, 0.5 - phi, 0.5]),
                   :H_1 => Vec3([phi, 0.5 + phi, 0.5]), :X => Vec3([0.0, 0.5, 0.5]),
                   :Y => Vec3([0.5, 0.0, 0.5]), :Z => Vec3([0.5, 0.5, 0.0]))
    path = [[:Γ, :Y, :C, :D, :X, :Γ, :Z, :D_1, :H, :C], [:C_1, :Z], [:X, :H_1], [:H, :Y],
            [:L, :Γ]]
    return (kpoints = kpoints, path = path)
end

function orcf3(a, b, c)
    zeta = (1 + a^2 / b^2 - a^2 / c^2) / 4
    eta = (1 + a^2 / b^2 + a^2 / c^2) / 4
    kpoints = Dict(:Γ => Vec3([0.0, 0.0, 0.0]), :A => Vec3([0.5, 0.5 + zeta, zeta]),
                   :A_1 => Vec3([0.5, 0.5 - zeta, 1 - zeta]), :L => Vec3([0.5, 0.5, 0.5]),
                   :T => Vec3([1, 0.5, 0.5]), :X => Vec3([0.0, eta, eta]),
                   :X_1 => Vec3([1, 1 - eta, 1 - eta]), :Y => Vec3([0.5, 0.0, 0.5]),
                   :Z => Vec3([0.5, 0.5, 0.0]))
    path = [[:Γ, :Y, :T, :Z, :Γ, :X, :A_1, :Y], [:X, :A, :Z], [:L, :Γ]]
    return (kpoints = kpoints, path = path)
end

function orci(a, b, c)
    zeta = (1 + a^2 / c^2) / 4
    eta = (1 + b^2 / c^2) / 4
    delta = (b^2 - a^2) / (4 * c^2)
    mu = (a^2 + b^2) / (4 * c^2)
    kpoints = Dict(:Γ => Vec3([0.0, 0.0, 0.0]), :L => Vec3([-mu, mu, 0.5 - delta]),
                   :L_1 => Vec3([mu, -mu, 0.5 + delta]),
                   :L_2 => Vec3([0.5 - delta, 0.5 + delta, -mu]),
                   :R => Vec3([0.0, 0.5, 0.0]), :S => Vec3([0.5, 0.0, 0.0]),
                   :T => Vec3([0.0, 0.0, 0.5]), :W => Vec3([0.25, 0.25, 0.25]),
                   :X => Vec3([-zeta, zeta, zeta]), :X_1 => Vec3([zeta, 1 - zeta, -zeta]),
                   :Y => Vec3([eta, -eta, eta]), :Y_1 => Vec3([1 - eta, eta, -eta]),
                   :Z => Vec3([0.5, 0.5, -0.5]))
    path = [[:Γ, :X, :L, :T, :W, :R, :X_1, :Z, :Γ, :Y, :S, :W], [:L_1, :Y], [:Y_1, :Z]]
    return (kpoints = kpoints, path = path)
end

function orcc(a, b, c)
    zeta = (1 + a^2 / b^2) / 4
    kpoints = Dict(:Γ => Vec3([0.0, 0.0, 0.0]), :A => Vec3([zeta, zeta, 0.5]),
                   :A_1 => Vec3([-zeta, 1 - zeta, 0.5]), :R => Vec3([0.0, 0.5, 0.5]),
                   :S => Vec3([0.0, 0.5, 0.0]), :T => Vec3([-0.5, 0.5, 0.5]),
                   :X => Vec3([zeta, zeta, 0.0]), :X_1 => Vec3([-zeta, 1 - zeta, 0.0]),
                   :Y => Vec3([-0.5, 0.5, 0]), :Z => Vec3([0.0, 0.0, 0.5]))
    path = [[:Γ, :X, :S, :R, :A, :Z, :Γ, :Y, :X_1, :A_1, :T, :Y], [:Z, :T]]
    return (kpoints = kpoints, path = path)
end

function hex()
    kpoints = Dict(:Γ => Vec3([0.0, 0.0, 0.0]), :A => Vec3([0.0, 0.0, 0.5]),
                   :H => Vec3([1.0 / 3.0, 1.0 / 3.0, 0.5]),
                   :K => Vec3([1.0 / 3.0, 1.0 / 3.0, 0.0]), :L => Vec3([0.5, 0.0, 0.5]),
                   :M => Vec3([0.5, 0.0, 0.0]))
    path = [[:Γ, :M, :K, :Γ, :A, :L, :H, :A], [:L, :M], [:K, :H]]
    return (kpoints = kpoints, path = path)
end

function rhl1(alpha)
    eta = (1 + 4 * cos(alpha)) / (2 + 4 * cos(alpha))
    nu = 3.0 / 4.0 - eta / 2.0
    kpoints = Dict(:Γ => Vec3([0.0, 0.0, 0.0]), :B => Vec3([eta, 0.5, 1.0 - eta]),
                   :B_1 => Vec3([1.0 / 2.0, 1.0 - eta, eta - 1.0]),
                   :F => Vec3([0.5, 0.5, 0.0]), :L => Vec3([0.5, 0.0, 0.0]),
                   :L_1 => Vec3([0.0, 0.0, -0.5]), :P => Vec3([eta, nu, nu]),
                   :P_1 => Vec3([1.0 - nu, 1.0 - nu, 1.0 - eta]),
                   :P_2 => Vec3([nu, nu, eta - 1.0]), :Q => Vec3([1.0 - nu, nu, 0.0]),
                   :X => Vec3([nu, 0.0, -nu]), :Z => Vec3([0.5, 0.5, 0.5]))
    path = [[:Γ, :L, :B_1], [:B, :Z, :Γ, :X], [:Q, :F, :P_1, :Z], [:L, :P]]
    return (kpoints = kpoints, path = path)
end

function rhl2(alpha)
    eta = 1 / (2 * tan(alpha / 2.0)^2)
    nu = 3.0 / 4.0 - eta / 2.0
    kpoints = Dict(:Γ => Vec3([0.0, 0.0, 0.0]), :F => Vec3([0.5, -0.5, 0.0]),
                   :L => Vec3([0.5, 0.0, 0.0]), :P => Vec3([1 - nu, -nu, 1 - nu]),
                   :P_1 => Vec3([nu, nu - 1.0, nu - 1.0]), :Q => Vec3([eta, eta, eta]),
                   :Q_1 => Vec3([1.0 - eta, -eta, -eta]), :Z => Vec3([0.5, -0.5, 0.5]))
    path = [[:Γ, :P, :Z, :Q, :Γ, :F, :P_1, :Q_1, :L, :Z]]
    return (kpoints = kpoints, path = path)
end

function mcl(b, c, beta)
    eta = (1 - b * cos(beta) / c) / (2 * sin(beta)^2)
    nu = 0.5 - eta * c * cos(beta) / b
    kpoints = Dict(:Γ => Vec3([0.0, 0.0, 0.0]), :A => Vec3([0.5, 0.5, 0.0]),
                   :C => Vec3([0.0, 0.5, 0.5]), :D => Vec3([0.5, 0.0, 0.5]),
                   :D_1 => Vec3([0.5, 0.5, -0.5]), :E => Vec3([0.5, 0.5, 0.5]),
                   :H => Vec3([0.0, eta, 1.0 - nu]), :H_1 => Vec3([0.0, 1.0 - eta, nu]),
                   :H_2 => Vec3([0.0, eta, -nu]), :M => Vec3([0.5, eta, 1.0 - nu]),
                   :M_1 => Vec3([0.5, 1 - eta, nu]), :M_2 => Vec3([0.5, 1 - eta, nu]),
                   :X => Vec3([0.0, 0.5, 0.0]), :Y => Vec3([0.0, 0.0, 0.5]),
                   :Y_1 => Vec3([0.0, 0.0, -0.5]), :Z => Vec3([0.5, 0.0, 0.0]))
    path = [[:Γ, :Y, :H, :C, :E, :M_1, :A, :X, :H_1], [:M, :D, :Z], [:Y, :D]]
    return (kpoints = kpoints, path = path)
end

function mclc1(a, b, c, alpha)
    zeta = (2 - b * cos(alpha) / c) / (4 * sin(alpha)^2)
    eta = 0.5 + 2 * zeta * c * cos(alpha) / b
    psi = 0.75 - a^2 / (4 * b^2 * sin(alpha)^2)
    phi = psi + (0.75 - psi) * b * cos(alpha) / c
    kpoints = Dict(:Γ => Vec3([0.0, 0.0, 0.0]), :N => Vec3([0.5, 0.0, 0.0]),
                   :N_1 => Vec3([0.0, -0.5, 0.0]),
                   :F => Vec3([1 - zeta, 1 - zeta, 1 - eta]),
                   :F_1 => Vec3([zeta, zeta, eta]), :F_2 => Vec3([-zeta, -zeta, 1 - eta]),
                   # :F_3 => Vec3([1 - zeta, -zeta, 1 - eta]),
                   :I => Vec3([phi, 1 - phi, 0.5]), :I_1 => Vec3([1 - phi, phi - 1, 0.5]),
                   :L => Vec3([0.5, 0.5, 0.5]), :M => Vec3([0.5, 0.0, 0.5]),
                   :X => Vec3([1 - psi, psi - 1, 0.0]), :X_1 => Vec3([psi, 1 - psi, 0.0]),
                   :X_2 => Vec3([psi - 1, -psi, 0.0]), :Y => Vec3([0.5, 0.5, 0.0]),
                   :Y_1 => Vec3([-0.5, -0.5, 0.0]), :Z => Vec3([0.0, 0.0, 0.5]))
    path = [[:Γ, :Y, :F, :L, :I], [:I_1, :Z, :F_1], [:Y, :X_1], [:X, :Γ, :N], [:M, :Γ]]
    return (kpoints = kpoints, path = path)
end

function mclc2(a, b, c, alpha)
    zeta = (2 - b * cos(alpha) / c) / (4 * sin(alpha)^2)
    eta = 0.5 + 2 * zeta * c * cos(alpha) / b
    psi = 0.75 - a^2 / (4 * b^2 * sin(alpha)^2)
    phi = psi + (0.75 - psi) * b * cos(alpha) / c
    kpoints = Dict(:Γ => Vec3([0.0, 0.0, 0.0]), :N => Vec3([0.5, 0.0, 0.0]),
                   :N_1 => Vec3([0.0, -0.5, 0.0]),
                   :F => Vec3([1 - zeta, 1 - zeta, 1 - eta]),
                   :F_1 => Vec3([zeta, zeta, eta]), :F_2 => Vec3([-zeta, -zeta, 1 - eta]),
                   :F_3 => Vec3([1 - zeta, -zeta, 1 - eta]),
                   :I => Vec3([phi, 1 - phi, 0.5]), :I_1 => Vec3([1 - phi, phi - 1, 0.5]),
                   :L => Vec3([0.5, 0.5, 0.5]), :M => Vec3([0.5, 0.0, 0.5]),
                   :X => Vec3([1 - psi, psi - 1, 0.0]), :X_1 => Vec3([psi, 1 - psi, 0.0]),
                   :X_2 => Vec3([psi - 1, -psi, 0.0]), :Y => Vec3([0.5, 0.5, 0.0]),
                   :Y_1 => Vec3([-0.5, -0.5, 0.0]), :Z => Vec3([0.0, 0.0, 0.5]))
    path = [[:Γ, :Y, :F, :L, :I], [:I_1, :Z, :F_1], [:N, :Γ, :M]]
    return (kpoints = kpoints, path = path)
end

function mclc3(a, b, c, alpha)
    mu = (1 + b^2 / a^2) / 4.0
    delta = b * c * cos(alpha) / (2 * a^2)
    zeta = mu - 0.25 + (1 - b * cos(alpha) / c) / (4 * sin(alpha)^2)
    eta = 0.5 + 2 * zeta * c * cos(alpha) / b
    phi = 1 + zeta - 2 * mu
    psi = eta - 2 * delta
    kpoints = Dict(:Γ => Vec3([0.0, 0.0, 0.0]), :F => Vec3([1 - phi, 1 - phi, 1 - psi]),
                   :F_1 => Vec3([phi, phi - 1, psi]),
                   :F_2 => Vec3([1 - phi, -phi, 1 - psi]), :H => Vec3([zeta, zeta, eta]),
                   :H_1 => Vec3([1 - zeta, -zeta, 1 - eta]),
                   :H_2 => Vec3([-zeta, -zeta, 1 - eta]), :I => Vec3([0.5, -0.5, 0.5]),
                   :M => Vec3([0.5, 0.0, 0.5]), :N => Vec3([0.5, 0.0, 0.0]),
                   :N_1 => Vec3([0.0, -0.5, 0.0]), :X => Vec3([0.5, -0.5, 0.0]),
                   :Y => Vec3([mu, mu, delta]), :Y_1 => Vec3([1 - mu, -mu, -delta]),
                   :Y_2 => Vec3([-mu, -mu, -delta]), :Y_3 => Vec3([mu, mu - 1, delta]),
                   :Z => Vec3([0.0, 0.0, 0.5]))
    path = [[:Γ, :Y, :F, :H, :Z, :I, :F_1], [:H_1, :Y_1, :X, :Γ, :N], [:M, :Γ]]
    return (kpoints = kpoints, path = path)
end

function mclc4(a, b, c, alpha)
    mu = (1 + b^2 / a^2) / 4.0
    delta = b * c * cos(alpha) / (2 * a^2)
    zeta = mu - 0.25 + (1 - b * cos(alpha) / c) / (4 * sin(alpha)^2)
    eta = 0.5 + 2 * zeta * c * cos(alpha) / b
    phi = 1 + zeta - 2 * mu
    psi = eta - 2 * delta
    kpoints = Dict(:Γ => Vec3([0.0, 0.0, 0.0]), :F => Vec3([1 - phi, 1 - phi, 1 - psi]),
                   :F_1 => Vec3([phi, phi - 1, psi]),
                   :F_2 => Vec3([1 - phi, -phi, 1 - psi]), :H => Vec3([zeta, zeta, eta]),
                   :H_1 => Vec3([1 - zeta, -zeta, 1 - eta]),
                   :H_2 => Vec3([-zeta, -zeta, 1 - eta]), :I => Vec3([0.5, -0.5, 0.5]),
                   :M => Vec3([0.5, 0.0, 0.5]), :N => Vec3([0.5, 0.0, 0.0]),
                   :N_1 => Vec3([0.0, -0.5, 0.0]), :X => Vec3([0.5, -0.5, 0.0]),
                   :Y => Vec3([mu, mu, delta]), :Y_1 => Vec3([1 - mu, -mu, -delta]),
                   :Y_2 => Vec3([-mu, -mu, -delta]), :Y_3 => Vec3([mu, mu - 1, delta]),
                   :Z => Vec3([0.0, 0.0, 0.5]))
    path = [[:Γ, :Y, :F, :H, :Z, :I], [:H_1, :Y_1, :X, :Γ, :N], [:M, :Γ]]
    return (kpoints = kpoints, path = path)
end

function mclc5(a, b, c, alpha)
    zeta = (b^2 / a^2 + (1 - b * cos(alpha) / c) / sin(alpha)^2) / 4
    eta = 0.5 + 2 * zeta * c * cos(alpha) / b
    mu = eta / 2 + b^2 / (4 * a^2) - b * c * cos(alpha) / (2 * a^2)
    nu = 2 * mu - zeta
    rho = 1 - zeta * a^2 / b^2
    omega = (4 * nu - 1 - b^2 * sin(alpha)^2 / a^2) * c / (2 * b * cos(alpha))
    delta = zeta * c * cos(alpha) / b + omega / 2 - 0.25
    kpoints = Dict(:Γ => Vec3([0.0, 0.0, 0.0]), :F => Vec3([nu, nu, omega]),
                   :F_1 => Vec3([1 - nu, 1 - nu, 1 - omega]),
                   :F_2 => Vec3([nu, nu - 1, omega]), :H => Vec3([zeta, zeta, eta]),
                   :H_1 => Vec3([1 - zeta, -zeta, 1 - eta]),
                   :H_2 => Vec3([-zeta, -zeta, 1 - eta]), :I => Vec3([rho, 1 - rho, 0.5]),
                   :I_1 => Vec3([1 - rho, rho - 1, 0.5]), :L => Vec3([0.5, 0.5, 0.5]),
                   :M => Vec3([0.5, 0.0, 0.5]), :N => Vec3([0.5, 0.0, 0.0]),
                   :N_1 => Vec3([0.0, -0.5, 0.0]), :X => Vec3([0.5, -0.5, 0.0]),
                   :Y => Vec3([mu, mu, delta]), :Y_1 => Vec3([1 - mu, -mu, -delta]),
                   :Y_2 => Vec3([-mu, -mu, -delta]), :Y_3 => Vec3([mu, mu - 1, delta]),
                   :Z => Vec3([0.0, 0.0, 0.5]))
    path = [[:Γ, :Y, :F, :L, :I], [:I_1, :Z, :H, :F_1], [:H_1, :Y_1, :X, :Γ, :N], [:M, :Γ]]
    return (kpoints = kpoints, path = path)
end

function tria()
    kpoints = Dict(:Γ => Vec3([0.0, 0.0, 0.0]), :L => Vec3([0.5, 0.5, 0.0]),
                   :M => Vec3([0.0, 0.5, 0.5]), :N => Vec3([0.5, 0.0, 0.5]),
                   :R => Vec3([0.5, 0.5, 0.5]), :X => Vec3([0.5, 0.0, 0.0]),
                   :Y => Vec3([0.0, 0.5, 0.0]), :Z => Vec3([0.0, 0.0, 0.5]))
    path = [[:X, :Γ, :Y], [:L, :Γ, :Z], [:N, :Γ, :M], [:R, :Γ]]
    return (kpoints = kpoints, path = path)
end

function trib()
    kpoints = Dict(:Γ => Vec3([0.0, 0.0, 0.0]), :L => Vec3([0.5, -0.5, 0.0]),
                   :M => Vec3([0.0, 0.0, 0.5]), :N => Vec3([-0.5, -0.5, 0.5]),
                   :R => Vec3([0.0, -0.5, 0.5]), :X => Vec3([0.0, -0.5, 0.0]),
                   :Y => Vec3([0.5, 0.0, 0.0]), :Z => Vec3([-0.5, 0.0, 0.5]))
    path = [[:X, :Γ, :Y], [:L, :Γ, :Z], [:N, :Γ, :M], [:R, :Γ]]
    return (kpoints = kpoints, path = path)
end
