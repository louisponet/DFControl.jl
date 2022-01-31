function pseudos(pseudoset, atsyms::Vector{Symbol}, fuzzy = ""; server="localhost")
    s = Servers.maybe_start(server)
    pseudo_paths = list_pseudoset(pseudoset, fuzzy, server=s)
    ps = Dict{Symbol,String}()
    for a in atsyms
        path = get(pseudo_paths, a, nothing)
        if path === nothing
            error("No pseudo for atom $a found in set $pseudoset.")
        end
        t = tempname()
        Servers.pull(s, path, t)
        ps[a] = read(t, String)
        rm(t)
    end
    return ps
end

function list_pseudoset(pseudoset, fuzzy=""; server="localhost")
    s = Servers.maybe_start(server)
    resp = HTTP.get(s, "/pseudos/$pseudoset", JSON3.write(fuzzy))
    if resp.status == 204
        error("No pseudoset $pseudoset found on Server $(s.name). Please first configure it using configure_pseudoset.")
    end
    return JSON3.read(resp.body, Dict{Symbol,String})
end

"""
    list_pseudosets(;server = "localhost")

Lists the pseudosets that have previously been set up.
"""
function list_pseudosets(;server = "localhost")
    s = Servers.maybe_start(server)
    return JSON3.read(HTTP.get(s, "/pseudos").body, Vector{String})
end

"""
    configure_pseudoset(set_name::String, dir::String; server = "localhost")

Reads the specified `dir` and sets up the pseudos for `set`.
"""
function configure_pseudoset(set_name::String, dir::String; server = "localhost")
    s = Servers.maybe_start(server)
    p = isabspath(dir) ? dir : joinpath(s, dir)
    n_pseudos = JSON3.read(HTTP.post(s, "/configure_pseudoset/" * p, 
                                     JSON3.write(set_name)).body, Int)
    @info "Configured $n_pseudos pseudos on Server $(s.name), found in dir $p."
end

"""
    rm_pseudoset!(set_name::String; server = "localhost")

Removes the pseudo set from the server.
"""
function rm_pseudoset!(set_name::String; server = "localhost")
    s = Servers.maybe_start(server)
    return HTTP.put(s, "/rm_pseudos", JSON3.write(set_name))
end
#---#

"""
    set_pseudos!(job::Job, set::String; server=job.server, specifier::String="", kwargs...)
    set_pseudos!(structure::Structure, set::String; server="localhost", specifier::String="", kwargs...)

Sets the pseudopotentials of the atoms inside the `structure` (or `job.structure`) to the ones of `set`.
`specifier` can be specified as a fuzzy match to select a specific pseudos if multiple pseudopotentials exist in the set.
Example:
```
set_pseudos!(job, "pbesol", specifier="rrkjus")
```
will select the pseudo file that contains "rrkjus" in the filename.

The pseudos will be searched for in the `server`.
"""
Structures.set_pseudos!(job::Job, args...; server=job.server,kwargs...) =
    Structures.set_pseudos!(job.structure, args...; server=server, kwargs...)
    
function Structures.set_pseudos!(str::Structures.Structure, set; server = "localhost", specifier = "", kwargs...)
    atsyms = unique(map(x -> x.element.symbol, str.atoms))
    return Structures.set_pseudos!(str, pseudos(set, atsyms, specifier, server=server);
                                   kwargs...)
end
    
