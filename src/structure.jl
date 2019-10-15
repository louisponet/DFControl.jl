
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
    addprojections!(projdict, str.atoms)
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
    create_supercell(structure::AbstractStructure, na::Int, nb::Int, nc::Int)

Takes a structure and creates a supercell from it with the given amount of extra cells (`na, nb, nc`) along the a, b, c direction.
"""
function create_supercell(structure::AbstractStructure, na::Int, nb::Int, nc::Int)
    orig_ats   = atoms(structure)
    orig_cell  = cell(structure)
    scale_mat  = diagm(0 => 1 .+ [na, nb, nc])
    new_cell   = orig_cell * scale_mat
    new_atoms  = eltype(orig_ats)[]
    for ia=0:na, ib=0:nb, ic=0:nc
        # if all((ia, ib, ic) .== 0)
        #     continue
        # end
        transl_vec = orig_cell * [ia, ib, ic]
        for at in orig_ats
	        cart_pos = position_cart(at) + transl_vec
	        cryst_pos = inv(new_cell) * cart_pos
            push!(new_atoms, Atom(at, cart_pos, Point3(cryst_pos)))
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

# #taken from Crystals.jl

# const DEFAULT_TOLERANCE = 1e-5

# """
#     potential_equivalents(cell::Mat3, tolerance::Real=$(DEFAULT_TOLERANCE))

# Figures out vectors in the lattice defined by the cell which have the same
# lengths as the unit vectors. These new vectors are
# potentially symmetrically equivalent.
# """
# function potential_equivalents(cell::Mat3{T}, tolerance::Real=DEFAULT_TOLERANCE) where {T}
#     V = volume(cell)
#     a0 = cell[:, 1]
#     a1 = cell[:, 2]
#     a2 = cell[:, 3]

#     lengths = reduce(+, cell .* cell, dims=1)
#     max_norm = mapreduce(i -> norm(cell[:, i]), max, 1:3, init=zero(T))

#     n0 = ceil(Int, max_norm * norm(cross(a1, a2)) / V)
#     n1 = ceil(Int, max_norm * norm(cross(a2, a0)) / V)
#     n2 = ceil(Int, max_norm * norm(cross(a0, a1)) / V)

#     gvectors = [Array{T, 1}[] for u in 1:length(lengths)]
#     for i in -n0:n0, j in -n1:n1, k in -n2:n2
#         g = cell * [i, j, k]
#         glength = g ⋅g
#         for (length, result) in zip(lengths, gvectors)
#             ustrip(abs(length - glength)) < tolerance && push!(result, g)
#         end
#     end

#     [hcat(gs...) for gs in gvectors]
# end

# """
#     point_group_operations(cell::Mat3)
#     point_group_operations(str::Structure)

# Determines the point group operations for a 3x3 unit cell.
# Rotations are determined from G-vector triplets with the same norm as the unit-cell vectors.

# Implementation taken from [ENUM](http://enum.sourceforge.net/).
# """
# function point_group(cell::Mat3{LengthType{T,V}}, tolerance::Real=DEFAULT_TOLERANCE) where {T, V}
#     avecs, bvecs, cvecs = potential_equivalents(cell, tolerance)
#     identity = Mat3(Matrix{T}(I, 3, 3))
#     result = [identity]

#     recip_cell = reciprocal(cell)
#     for i in 1:size(avecs, 2), j in 1:size(bvecs, 2), k in 1:size(cvecs, 2)

#         # (potential) rotation in cartesian coordinates
#         rotation = Mat3(hcat(avecs[:, i], bvecs[:, j], cvecs[:, k]))

#         # check operator is invertible
#         ustrip(volume(rotation)) ≥ tolerance || continue

#         # rotation in fractional coordinates
#         rotation_operator = rotation * recip_cell
#         any(ustrip.(abs.(rotation_operator .- identity)) .> tolerance) || continue

#         # check matrix is a rotation
#         if !all(ustrip.(abs.(rotation_operator * transpose(rotation_operator) .- identity)) .< tolerance)
#             continue
#         end

#         # check rotation not in list
#         index = findfirst(x -> all(ustrip.(abs.(x .- rotation_operator)) .< tolerance), result)
#         index ≠ nothing || push!(result, rotation_operator)
#     end
#     return result
# end

# # Taken from Crystals.jl Gruber
# function def_test(real::Real; tolerance::Real=default_tolerance)
#     real ≥ tolerance && return [0, 1]
#     real > -tolerance && return [1, 0]
#     return [0, 0]
# end

# function def_test(args; tolerance=default_tolerance)
#     result = [0, 0]
#     for u in args
#         result += def_test(u, tolerance=tolerance)
#     end
#     result
# end

# function def_gt_0(args...; tolerance=default_tolerance)
#     z₀, positive = def_test(args; tolerance=tolerance)
#     positive == 3 || (z₀ == 0 && positive == 1)
# end

# function n1_action(params::Vector, rinv::Matrix)
#     rinv[:, :] = rinv * [0 -1 0; -1 0 0; 0 0 -1]
#     params[1], params[2] = params[2], params[1]
#     params[4], params[5] = params[5], params[4]
# end

# function n2_action(params::Vector, rinv::Matrix)
#     rinv[:, :] = rinv * [-1 0 0; 0 0 -1; 0 -1 0]
#     params[2], params[3] = params[3], params[2]
#     params[5], params[6] = params[6], params[5]
# end

# function n3_action(params::Vector, rinv::Matrix; tolerance=DEFAULT_TOLERANCE)
#     i = params[4] ≤ -tolerance ? -1 : 1
#     j = params[5] ≤ -tolerance ? -1 : 1
#     k = params[6] ≤ -tolerance ? -1 : 1
#     rinv[:, :] = rinv * [i 0 0; 0 j 0; 0 0 k]
#     params[4:end] = abs.(params[4:end])
# end

# function n4_action(params::Vector, rinv::Matrix; tolerance=DEFAULT_TOLERANCE)
#     i = params[4] ≥ tolerance ? -1 : 1
#     j = params[5] ≥ tolerance ? -1 : 1
#     k = params[6] ≥ tolerance ? -1 : 1
#     update = diagm([i, j, k])
#     if i * j * k < 0
#         if k == 1 && params[6] > -tolerance
#             update[3, 3] = -1
#         elseif j == 1 && params[5] > -tolerance
#             update[3, 3] = -1
#         elseif i == 1 && params[4] > -tolerance
#             update[1, 1] = -1
#         elseif i * j * k == -1
#             error("Internal error")
#         end
#     end
#     rinv[:, :] = rinv * update
#     params[4:end] = -abs.(params[4:end])
# end

# function n5_action(params::Vector, rinv::Matrix; tolerance=DEFAULT_TOLERANCE)
#     pos_or_neg = params[4] > tolerance ? -1 : 1
#     rinv[:, :] = rinv * [1 0 0; 0 1 pos_or_neg; 0 0 1]
#     params[3] += params[2] + pos_or_neg * params[4]
#     params[4] += 2pos_or_neg * params[2]
#     params[5] += pos_or_neg * params[6]
# end

# function n6_action(params::Vector, rinv::Matrix; tolerance=DEFAULT_TOLERANCE)
#     pos_or_neg = params[5] > tolerance ? -1 : 1
#     rinv[:, :] = rinv * [1 0 pos_or_neg; 0 1 0; 0 0 1]
#     params[3] += params[1] + pos_or_neg * params[5]
#     params[4] += pos_or_neg * params[6]
#     params[5] += 2pos_or_neg * params[1]
# end

# function n7_action(params::Vector, rinv::Matrix; tolerance=DEFAULT_TOLERANCE)
#     pos_or_neg = params[6] > tolerance ? -1 : 1
#     rinv[:, :] = rinv * [1 pos_or_neg 0; 0 1 0; 0 0 1]
#     params[2] += params[1] + pos_or_neg * params[6]
#     params[4] += pos_or_neg * params[5]
#     params[6] += 2pos_or_neg * params[1]
# end

# function n8_action(params::Vector, rinv::Matrix)
#     rinv[:, :] = rinv * [1 0 1; 0 1 1; 0 0 1]
#     params[3] += sum(params[1:2]) + sum(params[4:end])
#     params[4] += 2params[2] + params[6]
#     params[5] += 2params[1] + params[6]
# end

# """
#     gruber(cell::Mat3;
#            tolerance::Real=DEFAULT_TOLERANCE, itermax::Int=50,
#            max_no_change::Int=10)

# Determines Gruber cell of an input cell.
# The Gruber cell is an optimal parameterization of a lattice, e.g. shortest
# cell-vectors and angles closest to 90 degrees. The new cell is in the same basis
# as the origin cell: no rotation has been incurred. The cell parameters are
# uniquely determined, even though the cell itself is not (certain symmetry
# operations leaving the structure unchanged may yield a more recognizable cell).
# If you want a unique Cartesian cell (in a different Cartesian basis), use
# the `niggly` algorithm.
# # Arguments
# * `cell::Mat3`: the input lattice cell-vectors. Cannot be singular.
# * `itermax::Int`: maximum number of iterations before bailing out
# * `tolerance::Real`: tolerance parameter when comparing real numbers
# * `max_no_change::Int`: Maximum number of times to go through algorithm
#   without changes before bailing out
# """
# function gruber(cell::Mat3{T};
#                              tolerance::Real=DEFAULT_TOLERANCE, itermax::Int=50,
#                              max_no_change::Int=10) where {T}
#     @assert volume(cell) > tolerance "Cell volume too small."

#     if itermax ≤ 0
#         itermax = typemax(itermax)
#     end
#     if max_no_change ≤ 0
#         max_no_change = typemax(max_no_change)
#     end

#     ε = tolerance
#     metric = transpose(cell) * cell
#     params = vcat(diag(metric), [2metric[2, 3], 2metric[1, 3], 2metric[1, 2]])
#     rinv = Matrix{T}(I, size(metric))
#     no_change, previous = 0, -params[1:3]
#     iteration = 0
#     for iteration in 1:itermax
#         condition0 =
#             (a, b, d, e) -> a ≥ b + ε || (abs(a - b) < ε && abs(d) ≥ abs(e) + ε)

#         condition0(params[1], params[2], params[4], params[5]) &&
#             n1_action(params, rinv)
#         condition0(params[2], params[3], params[5], params[6]) &&
#             (n2_action(params, rinv); continue)

#         if def_gt_0(params[4:end]...; tolerance=ε)
#             n3_action(params, rinv; tolerance=ε)
#         else
#             n4_action(params, rinv; tolerance=ε)
#             if all(abs.(previous .- params[1:3]) .< ε)
#                 no_change += 1
#             else
#                 no_change = 0
#             end
#             no_change < max_no_change || break
#         end

#         condition1 =
#             (d, b, e, f) -> abs(d) ≥ b + ε ||
#             (abs(d - b) < ε && 2e ≤ f - ε) ||
#             (abs(d + b) < ε && f ≤ -ε)
#         condition1(params[4], params[2], params[5], params[6]) &&
#             (n5_action(params, rinv; tolerance=ε); continue)
#         condition1(params[5], params[1], params[4], params[6]) &&
#             (n6_action(params, rinv; tolerance=ε); continue)
#         condition1(params[6], params[1], params[4], params[5]) &&
#             (n7_action(params, rinv; tolerance=ε); continue)

#         sum_no_c = sum(params[1:2]) + sum(params[4:end])
#         (
#             sum_no_c ≤ -ε ||
#             (abs(params[6] - params[1]) < ε && (2params[4] ≤ params[5] - ε)) ||
#             (abs(sum_no_c) < ε  && 2params[1] + 2params[5] + params[6] ≥ ε)
#         ) || break
#         n8_action(params, rinv)
#     end
#     iteration == itermax && error("Reached end of iteration without converging")
#     return Mat3(cell * rinv)
# end

# gruber(cell::Mat3{T}; kwargs...) where {T<:LengthType} =
#     gruber(ustrip.(cell); kwargs...) * unit(T)

# """
#     niggly(cell::AbstractMatrix; kwargs...)

# Determines a unique Cartesian cell equivalent to the input, with the shortest
# possible vectors and squarest angles. For an explanation of the parameters, see
# `gruber`. In practice, this function computes the cell-parameters of a `gruber` cell and
# then reconstructs the cell matrix. Hence, the result should be quite unique for any lattice
# representation, including any rotation of the underlying Cartesian coordinates.
# """
# niggly(cell::Mat3, kwargs...) = cell_parameters(cell_parameters(gruber(cell; kwargs...))...)


# """
#     cell_parameters(a, b, c, α=π/2, β=π/2, γₒ=π/2)

# Computes the cell matrix from the cell parameters [a, b, c, α, β, γ].
# """
# function cell_parameters(a::T, b, c, α=π/2, β=π/2, γₒ=π/2) where {T}
#     cx = cos(β)
#     cy = (cos(α) - cos(β)*cos(γₒ))/sin(γₒ)
#     z₀ = zero(T)
#     [a    b * cos(γₒ) c * cx;
#      z₀   b * sin(γₒ) c * cy;
#      z₀   z₀           c * √(1 - cx * cx - cy * cy)]
# end

# """
#     cell_parameters(cell::Mat3)
#     cell_parameters(s::Structure)
# Parameters (a, b, c, α, β, γ) of the input cell returned in a named tuple.
# """
# function cell_parameters(cell::Mat3)
#     G = transpose(cell) * cell
#     a, b, c = sqrt.(diag(G))
#     α = acos(0.5(G[2, 3] + G[3, 2])/(c * b))
#     β = acos(0.5(G[3, 1] + G[1, 3])/(a * c))
#     γₒ = acos(0.5(G[1, 2] + G[2, 1])/(a * b))
#     a, b, c, α, β, γₒ
# end
# cell_parameters(s::Structure) = cell_parameters(cell(s))

# function inner_translations_impl(cell::Mat3,
#                                  atoms::Vector{<:AbstractAtom};
#                                  tolerance::Real=DEFAULT_TOLERANCE)
#     grubcell = gruber(cell)
#     inv_grubcell = inv(grubcell)

#     # find species with minimum number of atoms
#     uats = unique(atoms)
#     minl, atom_index = findmin(length.(map(x -> filter(y -> isequal_species(x, y), atoms), uats)))

#     test_at = atoms[atom_index]
#     atom_center = position_cart(test_at)

#     translations = []
#     for at in atoms
#         isequal_species(at, test_at) && continue

#         translation = grubcell * into_voronoi(position_cart(at) - atom_center, grubcell)
#         all(abs.(translation) .< tolerance) && continue

#         is_mapping = true
#         for at1 in atoms
#             pos = position_cryst(at1) + translation
#             found = false
#             for at2 in atoms
#                 !isequal_species(at1, at2) && continue
#                 found = is_periodic(pos, position_cryst(at2), grubcell, tolerance)
#                 @show found
#                 found && break
#             end
#             (!found) && (is_mapping = false; break)
#         end
#         is_mapping && push!(translations, into_voronoi(translation, grubcell))
#     end
#     length(translations) == 0 && return similar(position_cryst.(atoms), (length(atoms), 0))
#     hcat(translations...)
# end

# is_periodic(pos1::Point3, pos2::Point3, cell::Mat3, tolerance::Real=DEFAULT_TOLERANCE) = 
#     all(abs.(mod.(pos1 .- pos2 .+ 0.5, 1) .- 0.5) .< tolerance)

# """
#     into_voronoi(position_cryst::Point3, cell::Mat3; extent::Integer=1)

# Folds position into first Brillouin zone of the input cell. Makes a well-meaning effort at
# returning the periodic image with the smallest possible norm. It recenter the atoms around
# the origin and then looks for the smallest periodic images within `-extent:extent` cells. If
# the cell is quite pathological, then the result will not be within the Voronoi cell.
# """
# function into_voronoi(position_cryst::Point3, cell::Mat3; extent::Integer=1)
#     zcentered = cell * origin_centered(position_cryst)
#     result = copy(zcentered)
#     n = norm(zcentered)
#     for i = -extent:extent, j = -extent:extent, k = -extent:extent
#         translation = cell * [i, j, k]
#         position = zcentered + translation
#         d = norm(position)
#         if d < n
#             result = position
#             n = d
#         end
#     end
#     return inv(cell) * result
# end

# """
#     origin_centered(position_cryst::Point3)

# Folds positions back to origin, such that each fractional component ``x_f`` is between
# ``-0.5\\leq x_f < 0.5``.
# """
# origin_centered(position_cryst::Point3) =
#     mod.(position_cryst .+ 0.5, -1) .+ 0.5

