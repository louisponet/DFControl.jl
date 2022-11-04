module Utils
    using Dates
    using RemoteHPC: Server
#---------------------Parse utils-----------------------------------------------------#
"""
Searches a directory for all files containing the key.
"""
searchdir(path::String, key) = joinpath.((path,), filter(x -> occursin(key, x), readdir(path)))
function searchdir(server::Server, path::String, key)
    if server.domain == "localhost"
        return searchdir(path, key)
    else
        return joinpath.((path,), filter(x -> occursin(key, x), readdir(server, path)))
    end
end
export searchdir

"""
Parse an array of strings into an array of a type.
"""
function parse_string_array(T::Type, array)
    return map(x -> (v = tryparse(T, x); v == nothing ? 0.0 : v), array)
end
export parse_string_array

"""
Parse a line for occurrences of type T.
"""
parse_line(T::Type, line::String) = parse_string_array(T, split(line))
export parse_line

"""
Splits a line using arguments, then strips spaces from the splits.
"""
strip_split(line, args...) = strip.(split(line, args...))
export strip_split

function replace_multiple(str, replacements::Pair{String,String}...)
    tstr = deepcopy(str)
    for r in replacements
        tstr = replace(tstr, r)
    end
    return tstr
end
export replace_multiple

function cut_after(line, c)
    t = findfirst(isequal(c), line)
    if t == nothing
        return line
    elseif t == 1
        return ""
    else
        return line[1:t-1]
    end
end
export cut_after

#--------------------------------------------------------------------------------------#
function fort2julia(f_type)
    f_type = lowercase(f_type)
    if f_type == "real"
        return Float32
    elseif f_type == "real(kind=dp)"
        return Float64
    elseif f_type == "complex(kind=dp)"
        return Complex{Float64}
    elseif occursin("character", f_type)
        return String
    elseif f_type == "string"
        return String
    elseif f_type == "integer"
        return Int
    elseif f_type == "logical"
        return Bool
    elseif occursin(".D", f_type)
        return replace(f_type, "D" => "e")
    else
        return Nothing
    end
end
export fort2julia

"""
It's like filter()[1].
"""
function getfirst(f::Function, A)
    for el in A
        if f(el)
            return el
        end
    end
end
export getfirst

"""
    parse_block(f, types...; to_strip=',')

Takes the specified types and parses each line into the types.
When it finds a line where it cannot match all the types, it stops and returns  the parsed values.
The split and strip keywords let the user specify how to first split the line, then strip the splits from the strip char.
"""
function parse_block(f, types...; to_strip = ',')
    output = []
    len_typ = length(types)
    while !eof(f)
        line = strip.(split(readline(f)), to_strip)
        len_lin = length(line)
        if isempty(line)
            continue
        end
        i, j = 1, 1
        tmp = []
        while i <= len_typ && j <= len_lin
            typ = types[i]
            l   = line[j]
            try
                t = Meta.parse(l)
                if typeof(t) == typ
                    push!(tmp, t)
                    i += 1
                    j += 1
                else
                    j += 1
                end
            catch
                j += 1
            end
        end
        if length(tmp) < length(types)
            return output
        end
        push!(output, Tuple{types...}(tmp))
    end
    return output
end
export parse_block

yesterday() = today() - Day(1)
lastweek()  = today() - Week(1)
lastmonth() = today() - Month(1)

username() = read(`whoami`, String)
export yesterday, lastweek, lastmonth, username

function loaded_modules_string()
    all = string.(values(Base.loaded_modules))
    valid = String[]
    ks = keys(Pkg.project().dependencies)
    for p in all
        if p ∈ ks
            push!(valid, p)
        end
    end
    return """using $(join(valid, ", "))\n"""
end
export loaded_modules_string

function string2cmd(str::AbstractString)
    Cmd(string.(split(str)))
end

end
