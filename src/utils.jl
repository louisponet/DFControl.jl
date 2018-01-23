"""
Searches a directory for all files containing the key.

Input: path::String,
key::String
"""
search_dir(path::String, key) = filter(x -> contains(x, key), readdir(path))

"""
Parse an array of strings into an array of a type.

Input:  T::Type,
array::Array{String,1}
Return: Array{T,1}
"""
parse_string_array(T::Type, array) = map(x -> (v = tryparse(T, x); isnull(v) ? 0.0 : get(v)), array)

"""
Parse a line for occurrences of type T.

Input:  T::Type,
line::String
Return: Array{T,1}
"""
parse_line(T::Type, line::String) = parse_string_array(T, split(line))

"""
Mutatatively applies the fermi level to all eigvals in the band. If fermi is a quantum espresso scf output file it will try to find it in there.

Input:  band::Band,
fermi::Union{String,AbstractFloat}
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
Same as above but for an array of bands. Is this even necessary?
"""
function apply_fermi_level!{T<:Band}(bands::Array{T}, fermi)
    for band in bands
        apply_fermi_level!(band,fermi)
    end
end

"""
Same as above but not mutatatively.
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

function apply_fermi_level(bands::Array{T}, fermi) where T <: Band
    out = similar(bands)
    for (i, band) in enumerate(bands)
        out[i] = apply_fermi_level(band, fermi)
    end
    return out
end

"""
Makes sure that a directory string ends with "/".

Input:  directory::String
Return: String
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

strip_split(line, args...) = strip.(split(line, args...))

#Incomplete for now only QE flags are returned
"""
    print_qe_flags(namelist_symbol::Symbol)

Prints the possible Quantum Espresso input flags and their type for a given input namelist.
"""
function print_qe_flags(namelist_symbol::Symbol)
    for block in QEControlFlags
        if block.name == namelist_symbol
            display(block)
        end
    end
end

"""
print_qe_namelists()

Prints all the possible Quantum Espresso input namelists.
"""
function print_qe_namelists()
    for block in QEControlFlags
        dfprintln(block.name)
    end
end

const assets_dir = joinpath(@__DIR__, "../assets/")

const conversions = OrderedDict{Symbol,Float64}(:bohr2ang => 0.529177)

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

function convert_atoms2symdict(T::Type{Union{Dict, OrderedDict}}, atoms::Array{<:AbstractString, 1}, U=Float64)
    at_dict = T{Symbol, Array{Point3D{U}, 1}}()
    for line in atoms
        atsym, x, y, z = parse.(split(line))
        if !haskey(at_dict, atsym)
            at_dict[atsym] = [Point3D{U}(x, y, z)]
        else
            push!(at_dict[atsym], Point3D{U}(x, y, z))
        end
    end
    return at_dict
end

function convert_atoms2symdict(T::Union{Type{Dict},Type{OrderedDict}}, atoms::V, U=Float64) where V<:Union{Dict, OrderedDict}
    at_dict = T{Symbol, Array{Point3D{U}, 1}}()
    for (atsym, at) in atoms
        if !haskey(at_dict, atsym)
            at_dict[atsym] = [Point3D{U}(at[1].x, at[1].y, at[1].z)]
        end
        for pos in at[2:end]
            push!(at_dict[atsym],Point3D{U}(pos.x, pos.y, pos.z))
        end
    end
    return at_dict
end

function convert_atoms2symdict(T::Union{Type{Dict},Type{OrderedDict}}, atoms::Array{Pair{Symbol, Array{<:AbstractFloat, 1}}, 1}, U=Float64)
    at_dict = T{Symbol, Array{Point3D{U}, 1}}()
    for (atsym, at) in atoms
        if !haskey(at_dict, atsym)
            at_dict[atsym] = [Point3D{U}(at[1]...)]
        end
        for pos in at[2:end]
            push!(at_dict[atsym],Point3D{U}(pos...))
        end
    end
    return at_dict
end

function firstval(f::Funtion, A)
    for el in A
        if f(el)
            return el
        end
    end
end