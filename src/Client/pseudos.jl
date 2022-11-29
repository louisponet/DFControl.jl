using ..Structures: Pseudo

Base.@kwdef struct PseudoSet <: RemoteHPC.Storable
    name::String
    pseudos::Dict{Symbol, Vector{Pseudo}} = Dict{Symbol, Vector{Pseudo}}()
end
PseudoSet(name::String; kwargs...) = PseudoSet(; name=name, kwargs...)

RemoteHPC.storage_directory(pseudos::PseudoSet) = "pseudos"

function set_server!(set::PseudoSet, s::Server)
    for ps in values(set.pseudos)
        for p in ps
            p.server = s.name
        end
    end
    return set
end
    
function RemoteHPC.load(server, s::PseudoSet)
    uri = RemoteHPC.storage_uri(s)
    if RemoteHPC.exists(server, s) # asking for a stored item
        set = JSON3.read(JSON3.read(HTTP.get(server, uri).body, String), PseudoSet)
        return set_server!(set, server)
    else
        res = HTTP.get(server, HTTP.URI(path="/storage/", query= Dict("path"=> splitdir(HTTP.queryparams(uri)["path"])[1])))
        if !isempty(res.body)
            return JSON3.read(res.body, Vector{String})
        else
            @warn "No PseudoSets found. Use `configure_pseudoset` first."
            return String[]
        end
    end
end

"""
    configure_pseudoset(set_name::String, dir::String)

Reads the specified `dir` and sets up the pseudos for `set`.
"""
function configure_pseudoset(server::Server, set_name::String, dir::String)
    files = readdir(server, dir)
    pseudos = Dict{Symbol,Vector{Pseudo}}([el.symbol => Pseudo[]
                                           for el in Structures.ELEMENTS])
    for pseudo_string in files
        element = Symbol(titlecase(String(split(split(replace(pseudo_string, "-" => "_"), ".")[1], "_")[1])))
        if haskey(pseudos, element)
            push!(pseudos[element], Pseudo("", joinpath(dir, pseudo_string), ""))
        end
    end

    set = PseudoSet(set_name, pseudos)
    save(server, set)
    set_server!(set, server)
    return set
end

Structures.set_pseudos!(job::Job, set::PseudoSet, args...) = 
    Structures.set_pseudos!(job.structure, set, args...)
    
function Structures.set_pseudos!(str::Structures.Structure, set::PseudoSet, args...)
    @assert !isempty(set.pseudos) "Error, no pseudos in pseudoset, load it first."
    return Structures.set_pseudos!(str, set.pseudos, args...)
end
    
