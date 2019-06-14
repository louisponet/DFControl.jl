
abstract type AbstractStructure{T} end

mutable struct Structure{T <: AbstractFloat, AA<:AbstractAtom{T}} <: AbstractStructure{T}
    name ::AbstractString
    cell ::Mat3{T}
    atoms::Vector{AA}
    data ::Dict{Symbol, Any}
end

Structure(name, cell::Mat3{T}, atoms::Vector{Atom{T}}) where T <: AbstractFloat = Structure{T, Atom{T}}(name, cell, atoms, Dict{Symbol, Any}())
Structure(str::AbstractStructure{T}, atoms::Vector{AT}) where {T <: AbstractFloat, AT<:AbstractAtom{T}} = Structure{T, AT}(name(str), cell(str), atoms, data(str))
Structure(cell::Matrix{T}, atoms::Vector{Atom{T}}) where T <: AbstractFloat = Structure{T, Atom{T}}("NoName", cell, atoms, Dict{Symbol, Any}())
Structure() = Structure("NoName", eye(3), Atom[], Dict{Symbol, Any}())
Structure(cif_file::String; name="NoName") = cif2structure(cif_file, structure_name = name)

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
                    if fname in [:name, :element, :position, :magnetization]
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

function setpseudos!(structure::AbstractStructure, set::Symbol, specifier::String=""; kwargs...)
    for (i, at) in enumerate(atoms(structure))
        pseudo = getdefault_pseudo(name(at), set, specifier=specifier)
        if pseudo == nothing
            @warn "Pseudo for $(name(at)) at index $i not found in set $set."
        else
            setpseudo!(at, pseudo; kwargs...)
        end
    end
end

function setpseudos!(structure::AbstractStructure, atname::Symbol, set::Symbol, specifier::String=""; kwargs...)
    for (i, at) in enumerate(atoms(structure, atname))
        pseudo = getdefault_pseudo(name(at), set, specifier=specifier)
        if pseudo == nothing
            @warn "Pseudo for $(name(at)) at index $i not found in set $set."
        else
            setpseudo!(at, pseudo; kwargs...)
        end
    end
end

function setpseudos!(structure::AbstractStructure, at_pseudos::Pair{Symbol, String}...; kwargs...)
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
    new_cell   = scale_mat * orig_cell
    new_atoms  = deepcopy(orig_ats)
    for ia=0:na, ib=0:nb, ic=0:nc
        if all((ia, ib, ic) .== 0)
            continue
        end
        transl_vec = orig_cell'*[ia, ib, ic]
        for at in orig_ats
            push!(new_atoms, Atom(at, position(at)+transl_vec))
        end
    end
    return Structure(name(structure), Mat3(new_cell), new_atoms, data(structure))
end
