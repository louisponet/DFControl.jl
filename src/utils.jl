"""
Searches a directory for all files containing the key.
"""
search_dir(path::String, key) = filter(x -> contains(x, key), readdir(path))

"""
Parse an array of strings into an array of a type.
"""
parse_string_array(T::Type, array) = map(x -> (v = tryparse(T, x); isnull(v) ? 0.0 : get(v)), array)

"""
Parse a line for occurrences of type T.
"""
parse_line(T::Type, line::String) = parse_string_array(T, split(line))

"""
Splits a line using arguments, then strips spaces from the splits.
"""
strip_split(line, args...) = strip.(split(line, args...))

"""
Mutatatively applies the fermi level to all eigvals in the band. If fermi is a quantum espresso scf output file it will try to find it in there.
"""
function apply_fermi_level!(band::Band, fermi::Union{String,AbstractFloat})
    if typeof(fermi) == String
        fermi = read_fermi_from_qe_file(fermi)
    end
    for i = 1:size(band.eigvals)[1]
        band.eigvals[i] -= fermi
    end
end

"""
Applies the fermi level to all eigvals in the band. If fermi is a quantum espresso scf output file it will try to find it in there.
"""
function apply_fermi_level(band::Band, fermi)
    T = typeof(band.eigvals[1])
    if typeof(fermi) == String
        fermi = read_fermi_from_qe_file(fermi)
    end
    out = deepcopy(band)
    for i1 = 1:size(band.eigvals)[1]
        out.eigvals[i1] = band.eigvals[i1] - T(fermi)
    end
    return out
end

"""
Makes sure that a directory string ends with "/".
"""
function form_directory(directory::String)
    if directory[end] != '/'
        return directory * "/"
    else
        return directory
    end
end

"""
    gen_k_grid(na, nb, nc, input, T=Float64)

Returns an array of k-grid points that are equally spaced, input can be either `:wan` or `:nscf`, the returned grids are appropriate as inputs for wannier90 or an nscf calculation respectively.
"""
function gen_k_grid(na, nb, nc, input, T=Float64)
    if input == :wan || typeof(input) == WannierInput
        return reshape([T[a, b, c] for a in collect(linspace(0, 1, na + 1))[1:end - 1], b in collect(linspace(0, 1, nb + 1))[1:end - 1], c in collect(linspace(0, 1, nc + 1))[1:end - 1]],(na * nb * nc))
    elseif input == :nscf || typeof(input) == QEInput
        return reshape([T[a, b, c, 1 / (na * nb * nc)] for a in collect(linspace(0, 1, na + 1))[1:end - 1], b in collect(linspace(0, 1, nb + 1))[1:end - 1], c in collect(linspace(0, 1, nc + 1))[1:end - 1]], (na * nb * nc))
    end
end

function fort2julia(f_type)
    f_type = lowercase(f_type)
    if f_type == "real"
        return Float32
    elseif f_type == "real(kind=dp)"
        return Float64
    elseif f_type == "complex(kind=dp)"
        return Complex{Float64}
    elseif contains(f_type, "character")
        return String
    elseif f_type == "string"
        return String
    elseif f_type == "integer"
        return Int
    elseif f_type == "logical"
        return Bool
    elseif contains(f_type,".D")
        return replace(f_type, "D", "e")
    end
end

function read_block(f, startstr::String, endstr::String)
    block = [startstr]
    line = readline(f)
    while !eof(f)
        while !contains(line, endstr)
            line = readline(f)
            push!(block, line)
        end
        return block
    end
    error("Block not found: start = $startstr, end = $endstr.")
end

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

function parse_block(f, types...; to_strip=',')
    output = []
    len_typ = length(types)
    while !eof(f)
        line = strip.(split(readline(f)), to_strip)
        len_lin = length(line)
        if isempty(line)
            continue
        end
        i,j = 1,1
        tmp = []
        while i <= len_typ && j <= len_lin
            typ = types[i]
            l   = line[j]
            try
                t   = parse(l)
                if typeof(t) == typ
                    push!(tmp, t)
                    i+=1
                    j+=1
                else
                    j+=1
                end
            catch
                j+=1
            end
        end
        if length(tmp) < length(types)
            return output
        end
        push!(output, Tuple{types...}(tmp))
    end
    return output
end
