
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

"""
sets the projections of the specified atoms.
"""
function setprojections!(str::Structure, projections...)
    projdict = Dict(projections)
    for at in unique_atoms(str.atoms)
        if !haskey(projdict, id(at))
            projdict[id(at)] = [proj.orb for proj in projections(at)]
        end
    end
    empty_projections!(str)
    addprojections!(projdict, str.atoms)
end

function empty_projections!(str::Structure)
    for at in str.atoms
        empty!(projections(at))
    end
end

"Takes a vector of structures and merges all the attributes of the atoms."
function merge_structures(structures::Vector{Union{<:AbstractStructure, Void}})
    nonvoid = filter(x -> x != nothing, structures)
    out_structure = nonvoid[1]
    for structure in nonvoid[2:end]
        for (at1, at2) in zip(out_structure.atoms, structure.atoms)
            for name in fieldnames(typeof(at1))
                if name in [:id, :element, :position, :pseudo]
                    continue
                end
                if !isdefined(at2, name)
                    continue
                else
                    setfield!(at1, name, getfield(at2,name))
                end
            end
        end
    end
    return out_structure
end

"Uses cif2cell to parse a cif file, then returns the parsed structure."
function cif2structure(cif_file::String; structure_name="NoName")
    tmpdir = tempdir()
    tmpfile = joinpath(tmpdir, "tmp.in")
    @assert splitext(cif_file)[2] == ".cif" error("Please specify a valid cif input file")
    run(`$cif2cellpath $cif_file -p quantum-espresso -o $tmpfile`)

    bla, structure = read_qe_input(tmpfile, structure_name = structure_name)
    rm(tmpfile)
    return structure
end
