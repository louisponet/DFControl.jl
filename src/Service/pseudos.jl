#We used this function to always keep the saved pseudos updated
function with_pseudos(f::Function)
    pp = config_path("pseudos.jld2")
    pseudos = ispath(pp) ? DFC.JLD2.load(pp)["pseudos"] : Dict{String, Dict{Symbol, Vector{Pseudo}}}()
    t = f(pseudos)
    DFC.JLD2.save(config_path("pseudos.jld2"), "pseudos", pseudos)
    return t
end

function pseudos(set::String, fuzzy::String)
    @info fuzzy
    with_pseudos() do all_sets
        out = Dict{Symbol, DFC.Pseudo}()
        pseudos = get(all_sets, set, nothing)
        pseudos === nothing && return nothing
        @info length(pseudos)
        for (k, ps) in pseudos
            p = getfirst(x->occursin(fuzzy, x.name), ps)
            if p !== nothing
                out[k] = p
            end
        end
        return out
    end
end

function rm_pseudos!(set::String)
    with_pseudos() do pseudos
        pop!(pseudos, set, nothing)
        return nothing
    end
end

"""
    configure_pseudos(set_name::String, dir::String)

Reads the specified `dir` and sets up the pseudos for `set`.
"""
function configure_pseudos(set_name::String, dir::String)
    with_pseudos() do pseudos
        files = readdir(dir)
        pseudos[set_name] = Dict{Symbol, Vector{Pseudo}}([el.symbol => Pseudo[] for el in DFC.ELEMENTS])
        for pseudo_string in files
            element = Symbol(titlecase(String(split(split(pseudo_string, ".")[1], "_")[1])))
            if haskey(pseudos[set_name], element)
                push!(pseudos[set_name][element], Pseudo(pseudo_string, splitdir(dir)[1]))
            end
        end
        return length(files)
    end
end

"""
    pseudo_sets()

Lists all pseudo_sets.
"""
pseudo_sets() = with_pseudos(x -> keys(x))