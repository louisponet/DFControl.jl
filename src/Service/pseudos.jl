#We used this function to always keep the saved pseudos updated
function with_pseudos(f::Function)
    pp = config_path("pseudos.jld2")
    pseudos = ispath(pp) ? JLD2.load(pp)["pseudos"] :
              Dict{String,Dict{Symbol,Vector{String}}}()
    t = f(pseudos)
    JLD2.jldsave(config_path("pseudos.jld2"); pseudos=pseudos)
    return t
end

function pseudos(set::String, fuzzy::String)
    with_pseudos() do all_sets
        out = Dict{Symbol,String}()
        pseudos = get(all_sets, set, nothing)
        pseudos === nothing && return nothing
        @info length(pseudos)
        for (k, ps) in pseudos
            p = getfirst(x -> occursin(fuzzy, x), ps)
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
    configure_pseudoset(set_name::String, dir::String)

Reads the specified `dir` and sets up the pseudos for `set`.
"""
function configure_pseudoset(set_name::String, dir::String)
    with_pseudos() do pseudos
        files = readdir(dir)
        pseudos[set_name] = Dict{Symbol,Vector{String}}([el.symbol => String[]
                                                         for el in Structures.ELEMENTS])
        for pseudo_string in files
            element = Symbol(titlecase(String(split(split(pseudo_string, ".")[1], "_")[1])))
            if haskey(pseudos[set_name], element)
                push!(pseudos[set_name][element], joinpath(dir, pseudo_string))
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
