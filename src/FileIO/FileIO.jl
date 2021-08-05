module FileIO
# This module handles all the io of files, and parsing.
using DelimitedFiles, Dates, LinearAlgebra, UnitfulAtomic
using UnitfulAtomic: ustrip, uconvert
using ..DFControl
using ..Utils
using ..Calculations
using ..Structures
using ..Jobs

include("qe.jl")
include("wannier.jl")
include("job.jl")

function parse_file(filename::AbstractString, parse_funcs::Vector{<:Pair{String}};
                    extra_parse_funcs::Vector{<:Pair} = Pair{String,Function}[])
    out = Dict{Symbol,Any}()
    open(filename, "r") do f
        lc = 0
        while !eof(f)
            line = strip(readline(f))
            lc += 1
            if isempty(line)
                continue
            end
            for pf in (parse_funcs, extra_parse_funcs)
                func = getfirst(x -> occursin(x[1], line), pf)
                func === nothing && continue
                try
                    func[2](out, line, f)
                catch
                    @warn "File corruption or parsing error detected executing parse function \n$(func[2]) in file $filename at line $lc: \"$line\".\nTrying to continue smoothly."
                end
            end
        end
    end
    return out
end

function cif2structure(cif_file::String; structure_name = "NoName")
    tmpdir = dirname(cif_file)
    tmpfile = joinpath(tmpdir, "tmp.in")
    @assert splitext(cif_file)[2] == ".cif" error("Please specify a valid cif calculation file")
    run(`$pythonpath $cif2cellpath $cif_file --no-reduce -p quantum-espresso -o $tmpfile`)

    bla, structure = qe_read_calculation(tmpfile; structure_name = structure_name)
    rm(tmpfile)
    return structure
end

function Base.parse(::Type{NamedTuple{names,types}},
                    spl::Vector{<:AbstractString}) where {names,types}
    @assert sum(length.(types.parameters)) == length(spl)
    tbuf = Vector{Union{types.parameters...}}(undef, length(types.parameters))
    counter = 1
    for (i, typ) in enumerate(types.parameters)
        if typ <: AbstractVector
            l = length(typ)
            tbuf[i] = parse(typ, spl[counter:counter+l-1])
            counter += l
        else
            tbuf[i] = parse(typ, spl[counter])
            counter += 1
        end
    end
    pstring = (tbuf...,)
    return NamedTuple{names}(pstring)
end
Base.parse(::Type{T}, s::AbstractString) where {T<:NamedTuple} = parse(T, split(s))

function Base.parse(::Type{Point{N,T}}, spl::Vector{<:AbstractString}) where {N,T}
    @assert N == length(spl)
    return Point{N,T}(parse.(T, spl))
end
Base.parse(::Type{T}, s::AbstractString) where {T<:Point} = parse(T, split(s))

end
