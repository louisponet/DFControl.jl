
abstract type AbstractStructure{T} end

mutable struct Structure{T <: AbstractFloat, AA<:AbstractAtom{T}} <: AbstractStructure{T}
    name ::AbstractString
    cell ::Mat3{T}
    atoms::Vector{AA}
    data ::Dict{Symbol, Any}
end

Structure(name, cell::Mat3{T}, atoms::Vector{Atom{T}}) where T <: AbstractFloat = Structure{T, Atom{T}}(name, cell, atoms, Dict{Symbol, Any}())
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
        id(at) == atsym && push!(out, at)
    end
    return out
end
atoms(str::AbstractStructure) = structure(str).atoms
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
        if !haskey(projdict, id(at))
            projdict[id(at)] = [proj.orb for proj in projections(at)]
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
            if at1==at2
                for name in fieldnames(typeof(at1))
                    if name in [:id, :element, :position]
                        continue
                    end
                    field =getfield(at2, name)
                    if field == nothing || isempty(field)
                        continue
                    else
                        setfield!(at1, name, getfield(at2,name))
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

function setpseudos!(structure::AbstractStructure, set, specifier=nothing)
    for (i, at) in enumerate(atoms(structure))
        pseudo = getdefault_pseudo(id(at), set, specifier=specifier)
        if pseudo == nothing
            @warn "Pseudo for $(id(at)) at index $i not found in set $set."
        else
            setpseudo!(at, pseudo)
        end
    end
end
function setpseudos!(structure::AbstractStructure, at_pseudos::Pair{Symbol, String}...)
    for (atsym, pseudo) in at_pseudos
        for at in atoms(structure, atsym)
            setpseudo!(at, pseudo)
        end
    end
end
