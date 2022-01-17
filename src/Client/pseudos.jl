function pseudos(server, pseudoset, atsyms::Vector{Symbol}, fuzzy = "")
    s = Servers.maybe_start_server(server)
    pseudo_paths = list_pseudoset(s, pseudoset, fuzzy)
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

function list_pseudoset(server, pseudoset, fuzzy="")
    s = Servers.maybe_start_server(server)
    resp = HTTP.get(s, "/pseudos/$pseudoset", [], JSON3.write(fuzzy))
    if resp.status == 204
        error("No pseudoset $pseudoset found on Server $(s.name). Please first configure it using configure_pseudoset.")
    end
    return JSON3.read(resp.body, Dict{Symbol,String})
end

"""
    list_pseudosets(server = "localhost")

Lists the pseudosets that have previously been set up.
"""
function list_pseudosets(server = "localhost")
    s = Servers.maybe_start_server(server)
    return JSON3.read(HTTP.get(s, "/pseudos").body, Vector{String})
end

"""
    configure_pseudoset(set_name::String, dir::String; server = "localhost")

Reads the specified `dir` and sets up the pseudos for `set`.
"""
function configure_pseudoset(set_name::String, dir::String; server = "localhost")
    s = Servers.maybe_start_server(server)
    p = isabspath(dir) ? dir : joinpath(s, dir)
    n_pseudos = JSON3.read(HTTP.post(s, "/configure_pseudoset/" * p, [],
                                     JSON3.write(set_name)).body, Int)
    @info "Configured $n_pseudos pseudos on Server $(s.name), found in dir $p."
end

"""
    rm_pseudoset!(set_name::String; server = "localhost")

Removes the pseudo set from the server.
"""
function rm_pseudoset!(set_name::String; server = "localhost")
    s = Servers.maybe_start_server(server)
    return HTTP.put(s, "/rm_pseudos", [], JSON3.write(set_name))
end
#---#

"""
    set_pseudos!(job::Job, set::Symbol, specifier::String=""; kwargs...)
    set_pseudos!(structure::Structure, set::Symbol, specifier::String=""; kwargs...)
    set_pseudos!(job::Job, atsym::Symbol, set::Symbol, specifier::String=""; kwargs...)
    set_pseudos!(structure::Structure, atsym::Symbol, set::Symbol, specifier::String=""; kwargs...)

Sets the pseudopotentials of the atoms inside the `structure` (or `job.structure`) to the ones of `set`.
`specifier` can be specified to select a specific pseudo if multiple pseudopotentials
for a given element exist in the set.
Example:
```
set_pseudos!(job, :pbesol, "rrkjus")
```
will select the pseudo file that contains "rrkjus" in the filename.

If `atsym` is used, only the pseudos of the atoms with that name will be set.

    set_pseudos!(job::Job, at_pseudos::Pair{Symbol, Pseudo}...; kwargs...)
    set_pseudos!(structure::Structure, at_pseudos::Pair{Symbol, Pseudo}...; kwargs...)

Convenience function that allows to set pseudopotentials for multiple atom types at the same time.
e.g. `set_pseudos!(job, :Si => getdefault_pseudo(:Si, :sssp)
"""
function Structures.set_pseudos!(job::Job, set; server=job.server, specifier::String = "", kwargs...)
    atsyms = unique(map(x -> x.element.symbol, job.structure.atoms))
    return Structures.set_pseudos!(job.structure, pseudos(server, set, atsyms, specifier);
                                   kwargs...)
end
