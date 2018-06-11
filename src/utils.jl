"""
Searches a directory for all files containing the key.
"""
searchdir(path::String, key) = filter(x -> contains(x, key), readdir(path))

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
        fermi = read_fermi_from_qe_output(fermi)
    end
    for i = 1:size(band.eigvals)[1]
        band.eigvals[i] -= fermi
    end
end

function Emax(Emin, nbnd, bands)
    nbndfound = 0
    max = 0
    for b in bands
        if minimum(b.eigvals) >= Emin && nbndfound <= nbnd
            nbndfound += 1
            #maximum of allowed frozen window is the minimum of the first band>nbnd
            max = minimum(b.eigvals)-0.005
        end
    end

    nbndfound <= nbnd && error("num_bands ($nbnd) starting from Emin=$Emin exceeds the number of bands ($nbndfound).")
    return max
end

function wanenergyranges(Emin, nbnd, bands, Epad=5)
    max = Emax(Emin, nbnd, bands)
    (Emin - Epad, Emin, max, max + Epad)
end

"""
Applies the fermi level to all eigvals in the band. If fermi is a quantum espresso scf output file it will try to find it in there.
"""
function apply_fermi_level(band::Band, fermi)
    T = typeof(band.eigvals[1])
    if typeof(fermi) == String
        fermi = read_fermi_from_qe_output(fermi)
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
function formdirectory(directory::String)
    if directory[end] != '/'
        return directory * "/"
    else
        return directory
    end
end

"""
    kgrid(na, nb, nc, input)

Returns an array of k-grid points that are equally spaced, input can be either `:wan` or `:nscf`, the returned grids are appropriate as inputs for wannier90 or an nscf calculation respectively.
"""
kgrid(na, nb, nc, ::Type{QE}) = reshape([[a, b, c, 1 / (na * nb * nc)] for a in collect(linspace(0, 1, na + 1))[1:end - 1], b in collect(linspace(0, 1, nb + 1))[1:end - 1], c in collect(linspace(0, 1, nc + 1))[1:end - 1]], (na * nb * nc))
kgrid(na, nb, nc, ::Type{Wannier90}) = reshape([[a, b, c] for a in collect(linspace(0, 1, na + 1))[1:end - 1], b in collect(linspace(0, 1, nb + 1))[1:end - 1], c in collect(linspace(0, 1, nc + 1))[1:end - 1]],(na * nb * nc))
kgrid(na, nb, nc, input::Symbol) = input==:wan ? kgrid(na, nb, nc, Wannier90) : kgrid(na, nb, nc, QE)
kgrid(na, nb, nc, input::DFInput{T}) where T = kgrid(na, nb, nc, T)

kakbkc(kgrid) = length.(unique.([[n[i] for n in kgrid] for i=1:3]))

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
    else
        return Void
    end
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

"""
    parse_block(f, types...; to_strip=',')

Takes the specified types and parses each line into the types.
When it finds a line where it cannot match all the types, it stops and returns  the parsed values.
The split and strip keywords let the user specify how to first split the line, then strip the splits from the strip char.
"""
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

convertflag(flagtype, flagvalue) = length(flagvalue) > 1 ? convert.(flagtype, flagvalue) : convert(flagtype, flagvalue)


macro undoable(func)
    if isa(func, Expr) && func.head === :function

        if func.args[1].args[2].head !== parameters
            # println(func.args[1].args[2])
            @assert func.args[1].args[2].args[2] == :DFJob error("This macro only works for `DFJob`.")
            j = func.args[1].args[2].args[1]
            bod = quote UNDO_JOBS[$(j).id]=deepcopy($(j))
                $(func.args[2])
                end
            Core.@__doc__ ret = Expr(func.head, func.args[1], bod)
            return ret
        else
            @assert func.args[1].args[3].args[2] == :DFJob error("This macro only works for `DFJob`.")
            j = func.args[1].args[3].args[1]
            bod = quote UNDO_JOBS[$(j).id]=deepcopy($(j))
                $(func.args[2])
                end
            Core.@__doc__ ret = Expr(func.head, func.args[1], bod)
            return ret
        end
    end
end
