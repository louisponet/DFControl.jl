const PSEUDO_DIR = DFC.config_path("pseudos")

function pseudos(set::String, fuzzy::String)
    if !(set in pseudo_sets())
        error("$set not found in pseudosets.")
    end
    out = Dict{Symbol,String}()
    
    pseudos = JSON3.read(read(joinpath(PSEUDO_DIR, set * ".json"), String), Dict{Symbol, Vector{String}})
    for (k, ps) in pseudos
        p = getfirst(x -> occursin(fuzzy, x), ps)
        if p !== nothing
            out[k] = p
        end
    end
    return out
end

function rm_pseudos!(set::String)
    p = joinpath(PSEUDO_DIR, set * ".json")
    if ispath(p)
        rm(p)
    end
end

"""
    configure_pseudoset(set_name::String, dir::String)

Reads the specified `dir` and sets up the pseudos for `set`.
"""
function configure_pseudoset(set_name::String, dir::String)
    files = readdir(dir)
    pseudos = Dict{Symbol,Vector{String}}([el.symbol => String[]
                                           for el in Structures.ELEMENTS])
    for pseudo_string in files
        element = Symbol(titlecase(String(split(split(pseudo_string, ".")[1], "_")[1])))
        if haskey(pseudos, element)
            push!(pseudos[element], joinpath(dir, pseudo_string))
        end
    end
    JSON3.write(joinpath(PSEUDO_DIR, set_name * ".json"), pseudos) 
    return length(pseudos)
end

"""
    pseudo_sets()

Lists all pseudo_sets.
"""
pseudo_sets() = map(x -> splitext(x)[1], readdir(PSEUDO_DIR))
