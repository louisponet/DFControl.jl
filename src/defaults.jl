"File with all the user defaults inside it"
# const default_file = abspath(homedir(),".julia","config","DFControl", "user_defaults.jl")
const default_file = occursin("cache", first(Base.DEPOT_PATH)) ? abspath(Base.DEPOT_PATH[2], "config","DFControl", "user_defaults.jl") : abspath(Base.DEPOT_PATH[1], "config","DFControl", "user_defaults.jl")
const default_pseudodirs = Dict{Symbol, String}()
const default_pseudos = Dict{Symbol, Dict{Symbol, Vector{Pseudo}}}()
const default_jobheader = [""]

for el in ELEMENTS
    default_pseudos[el.symbol] = Dict{Symbol, Vector{Pseudo}}()
end

const default_server = "localhost"

init_defaults(filename::String) = include(default_file)

"""
    setdefault_pseudodir(pseudo_symbol::Symbol, dir::String)

Adds an entry inside the `default_pseudodirs` with flag `pseudo_symbol`, and adds it to the `user_defaults.jl` file.
"""
function setdefault_pseudodir(pseudo_symbol::Symbol, dir::String)
    default_pseudodirs[pseudo_symbol] = dir
    expr2file(default_file, :(default_pseudodirs[$(QuoteNode(pseudo_symbol))] = $dir))
end

"""
    removedefault_pseudodir(pseudo_symbol::Symbol)

Removes entry with flag `pseudo_symbol` from the `default_pseudodirs` and `user_defaults.jl` file.
"""
function removedefault_pseudodir(pseudo_symbol::Symbol)
    if haskey(default_pseudodirs, pseudo_symbol)
        delete!(default_pseudodirs, pseudo_symbol)
        rm_expr_lhs(default_file, :(default_pseudodirs[$(QuoteNode(pseudo_symbol))]))
        # removedefault_pseudos(pseudo_symbol)
    end
end

"""
    removedefault_pseudos(pseudo_symbol::Symbol)

Removes all pseudo entries with flag `pseudo_symbol` from the `default_pseudos`.
"""
function removedefault_pseudos(pseudo_symbol::Symbol)
    found = false
    for (at, pseudos) in default_pseudos
        if haskey(pseudos, pseudo_symbol)
            delete!(pseudos, pseudo_symbol)
            rm_expr_lhs(default_file, :(default_pseudos[$(QuoteNode(at))][$(QuoteNode(pseudo_symbol))]))
            found = true
        end
    end
    if found
        removedefault_pseudodir(pseudo_symbol)
    end
end

"""
    setdefault_server(server::String)

Sets the default server variable, and also adds it to the `user_defaults.jl` file.
"""
function setdefault_server(server::String)
    default_server = server
    expr2file(default_file, :(default_server = $server))
end

"""
    getdefault_server()

Returns the default server if it's defined. If it is not defined return "".
"""
getdefault_server() = DFControl.default_server

"""
    getdefault_pseudodirs()

Returns the default pseudo dirs if it's defined. If it is not defined return nothing.
"""
getdefault_pseudodirs() = DFControl.default_pseudodirs

getdefault_pseudodir(pseudoset) = haskey(getdefault_pseudodirs(), pseudoset) ? getdefault_pseudodirs()[pseudoset] : nothing

"""
    configuredefault_pseudos(server = getdefault_server(), pseudo_dirs=getdefault_pseudodirs())

Reads the specified `default_pseudo_dirs` on the `default_server` and sets up the `default_pseudos` variable, and also adds all the entries to the `user_defaults.jl` file.
"""
function configuredefault_pseudos(;server = getdefault_server(), pseudo_dirs=getdefault_pseudodirs())
    if server == ""
        error("Either supply a valid server string or setup a default server through 'setdefault_server!()'.")
    end

    if pseudo_dirs == nothing
        error("Either supply valid pseudo directories or setup a default pseudo dir through 'setdefault_pseudodir()'.")
    end

    outputs = Dict{Symbol, Vector{String}}()
    for (name, directory) in pseudo_dirs
        outputs[name] = server == "localhost" ? readdir(directory) : split(read(`ssh -t $server ls $directory`, String), "\n")
    end

    elsyms = Symbol[el.symbol for el in ELEMENTS]

    for (name, pseudo_string) in outputs
        pseudos = filter(x -> x != "", pseudo_string)
        i = 1
        while i <= length(pseudos)
            pseudo  = pseudos[i]
            element = Symbol(titlecase(String(split(split(pseudo, ".")[1], "_")[1])))
            if element in elsyms
                t  = Pseudo[Pseudo(pseudo, pseudo_dirs[name])]
                j = 1
                while j + i <= length(pseudos) && Symbol(split(pseudos[i + j],".")[1]) == element
                    push!(t, Pseudo(pseudos[i + j], pseudo_dirs[name]))
                    j += 1
                end
                i += j
                expr2file(default_file, :(default_pseudos[$(QuoteNode(element))][$(QuoteNode(name))] = $t))
                default_pseudos[element][name] = t
            else
                i+=1
            end
        end
    end
end

"""
    getdefault_pseudo(atom::Symbol, set=:default; specifier=nothing)

Returns the pseudo potential string linked to the atom.
"""
function getdefault_pseudo(atom::Symbol, set=:default; specifier="")
    if tryparse(Int, String(atom)[end:end]) != nothing
        pp_atom = Symbol(String(atom)[1:end-1])
    else
        pp_atom = atom
    end
    if haskey(default_pseudos[pp_atom], set)
        if specifier != ""
            return deepcopy(getfirst(x -> occursin(specifier, x.name), default_pseudos[pp_atom][set]))
        else
            return deepcopy(default_pseudos[pp_atom][set][1])
        end
    end
end

"""
    setdefault_jobheader(lines)

Sets the header that will get added to each job.tt file, if no other header was specified.
"""
function setdefault_jobheader(lines)
    expr = :(default_jobheader = $lines)
    expr2file(default_file,expr)
    default_jobheader = lines
end

function getdefault_jobheader()
    if isdefined(DFControl, :default_job_header)
        return default_job_header
    else
        return [""]
    end
end

function findspecifier(str, strs::Vector{<:AbstractString})
    tmp = Char[]
    i = 1
    for s in strs
        if s == str
            continue
        end
        for (ch1, ch2) in zip(str, s)
            if ch1 != ch2
                push!(tmp, ch1)
            elseif !isempty(tmp)
                break
            end
        end
        !isempty(tmp) && break
        i += 1
    end
    testout = join(tmp)
    for s in strs[i+1:end]
        if occursin(testout, s)
            return findspecifier(str, strs[i+1:end])
        end
    end
    return testout
end

function getpseudoset(elsym::Symbol, ps::Pseudo)
	str = ps.name
    for (key, val) in default_pseudos[elsym]
        if length(val) == 1
            str == val[1] && return key, ""
        else
            for v in val
                if v == str
                    return key, findspecifier(str, val)
                end
            end
        end
    end
    return :none, ""
end
