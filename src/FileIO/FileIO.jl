module FileIO
# This module handles all the io of files, and parsing.
using DelimitedFiles, Dates, LinearAlgebra, UnitfulAtomic, CodeTracking
using UnitfulAtomic.Unitful: ustrip, uconvert
using UnitfulAtomic: bohr
using RemoteHPC: Exec, exec
using ..Utils
using ..Calculations
using ..Calculations: Calculation
using ..Structures
using ..Structures: Pseudo
using ..Jobs
using ..DFControl: Point,Point3, Vec3, SVector, Mat3, Mat4, Band, TimingData

include("qe.jl")
include("wannier.jl")
include("julia.jl")

function parse_file(f::IO, parse_funcs::Vector{<:Pair{String}};out = Dict{Symbol,Any}(),
                    extra_parse_funcs::Vector{<:Pair} = Pair{String,Function}[])
    
    lc = 0
    while !eof(f) && !haskey(out, :finished)
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
                @warn "File corruption or parsing error detected executing parse function \n$(func[2]) in file $f at line $lc: \"$line\".\nTrying to continue smoothly."
            end
        end
    end
    return out
end

function parse_file(f::AbstractString, args...; kwargs...)
    open(f, "r") do file
        parse_file(file, args...;kwargs...)
    end
end

function cif2structure(cif_file::String; structure_name = "NoName")
    tmpdir = dirname(cif_file)
    tmpfile = joinpath(tmpdir, "tmp.in")
    @assert splitext(cif_file)[2] == ".cif" error("Please specify a valid cif calculation file")
    run(`$pythonpath $cif2cellpath $cif_file --no-reduce -p quantum-espresso -o $tmpfile`)

    bla, structure = qe_parse_calculation(tmpfile; structure_name = structure_name)
    rm(tmpfile)
    return structure
end

Base.length(::Type{<:Real}) = 1

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

function write_xsf(filename::AbstractString, structure::Structure)
    open(filename, "w") do f
        write(f, "CRYSTAL\n")
        c = ustrip.(structure.cell')
        write(f, "PRIMVEC\n")
        write(f, "$(c[1,1]) $(c[1,2]) $(c[1,3])\n")
        write(f, "$(c[2,1]) $(c[2,2]) $(c[2,3])\n")
        write(f, "$(c[3,1]) $(c[3,2]) $(c[3,3])\n")
        write(f, "CONVVEC\n")
        write(f, "$(c[1,1]) $(c[1,2]) $(c[1,3])\n")
        write(f, "$(c[2,1]) $(c[2,2]) $(c[2,3])\n")
        write(f, "$(c[3,1]) $(c[3,2]) $(c[3,3])\n")
        write(f, "PRIMCOORD\n")
        write(f, "$(length(structure.atoms)) 1\n")
        for at in structure.atoms
            n = at.element.symbol
            p = ustrip.(at.position_cart)
            write(f, "$n $(p[1]) $(p[2]) $(p[3])\n")
        end
    end
end

"""
    outputdata(calculation::Calculation, file; extra_parse_funcs=[], print=true, overwrite=true)

If `file` exists this will parse it and return a `Dict` with the parsed data.
If `overwrite=false` and `calculation.outputdata` is not empty, this will be returned instead of reparsing the
output file.

`extra_parse_funcs` should be `Vector{Pair{String,Function}}`, where the string will be used to match `occursin(str, line)`
for each line of the output file. If this returns `true`, the function will be called as `func(results_dict, line, file)`.
The purpose is to allow for additional parsing that is not implemented, or for temporary non-standard values that are printed
while working on the DFT code, e.g. for debugging.

!!! note
    This only works for files where not all information is pulled out for yet,
    e.g. projwfc.x outputs are fully parsed already.

Example (from src/qe/fileio.jl):
```julia

function qe_parse_nat(results, line, f)
    results[:nat] = parse(Int, split(line)[end])
end

outputdata(job["scf"],joinpath(job, job["scf"].outfile), extra_parse_funcs = ["number of atoms/cell" => qe_parse_nat])
```
"""
function outputdata(calculation::Calculation, files...;
                    extra_parse_funcs::Vector{<:Pair{String}} = Pair{String}[],
                    print = true, overwrite = true)
    t = readoutput(calculation, files...; parse_funcs = extra_parse_funcs)
    return t === nothing ?
                          parse_file(files[1],
                                            extra_parse_funcs) : t
end

function calculationparser(e::Exec)
    if Calculations.is_qe_exec(e)
        qe_parse_calculation
    elseif Calculations.is_wannier_exec(e)
        wan_parse_calculation
    elseif Calculations.is_julia_exec(e)
        julia_parse_calculation
    end
end

end
