name(data::InputData) = data.name

name(c::DFCalculation)          = c.name
dir(c::DFCalculation)           = c.dir
flags(c::DFCalculation)         = c.flags
set_dir!(c::DFCalculation, dir) = (c.dir = dir)
name_ext(c::DFCalculation, ext) = name(c) * ext
infilename(c::DFCalculation)    = c.infile
outfilename(c::DFCalculation)   = c.outfile
inpath(c::DFCalculation)        = joinpath(c, infilename(c))
outpath(c::DFCalculation)       = joinpath(c, outfilename(c))

isbands(calculation::DFCalculation)    = false
isnscf(calculation::DFCalculation)     = false
isscf(calculation::DFCalculation)      = false
isvcrelax(calculation::DFCalculation)  = false
isprojwfc(calculation::DFCalculation)  = false
ismagnetic(calculation::DFCalculation) = false
issoc(calculation::DFCalculation)      = false

searchdir(i::DFCalculation, glob) = joinpath.((i,), searchdir(dir(i), glob))

"""
    joinpath(calc::DFCalculation, path...) = joinpath(dir(calc), path...)
"""
Base.joinpath(c::DFCalculation, path...) = joinpath(dir(c), path...)
hasflag(c::DFCalculation, s::Symbol) = haskey(flags(c), s)
hasflag(c::DFCalculation, s) = false

function flag(calculation::DFCalculation, flag::Symbol)
    if hasflag(calculation, flag)
        return calculation.flags[flag]
    end
end

function Base.:(==)(d1::InputData, d2::InputData)
    return all(x -> getfield(d1, x) == getfield(d2, x), fieldnames(InputData))
end

Base.get(c::DFCalculation, args...) = get(flags(c), args...)

package(::DFCalculation{T}) where {T} = T

data(calculation::DFCalculation) = calculation.data

hasexec(calculation::DFCalculation, ex::AbstractString) = exec(calculation, ex) != nothing
set_flow!(calculation::DFCalculation, run) = calculation.run = run


#TODO review this!
outdata(calculation::DFCalculation) = calculation.outdata
function hasoutput(calculation::DFCalculation)
    return !isempty(outdata(calculation)) || hasoutfile(calculation)
end

hasoutfile(calculation::DFCalculation) = ispath(outpath(calculation))
hasnewout(calculation::DFCalculation, time) = mtime(outpath(calculation)) > time

function set_cutoffs!(calculation::DFCalculation, args...)
    @warn "Setting cutoffs is not implemented for package $(package(calculation))"
end

function hasoutput_assert(i::DFCalculation)
    @assert hasoutfile(i) "Please specify an calculation that has an outputfile."
end

function hasexec_assert(i::DFCalculation, exec::String)
    @assert hasexec(i, exec) "Please specify an calculation with $exec as it's executable."
end

function Base.:(==)(i1::DFCalculation, i2::DFCalculation)
    return all(x -> x in (:outdata, :dir, :run) ? true : getfield(i1, x) == getfield(i2, x),
               fieldnames(DFCalculation))
end

"""
    kgrid(na, nb, nc, calculation)

Returns an array of k-grid points that are equally spaced, calculation can be either `:wan` or `:nscf`, the returned grids are appropriate as calculations for wannier90 or an nscf calculation respectively.
"""
kgrid(na, nb, nc, ::DFCalculation{T}) where {T} = kgrid(na, nb, nc, T)

include("qe/calculation.jl")
include("wannier90/calculation.jl")
include("elk/calculation.jl")
