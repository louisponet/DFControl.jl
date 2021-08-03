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

# Interface Functions
isbands(calculation::DFCalculation)    = false
isnscf(calculation::DFCalculation)     = false
isscf(calculation::DFCalculation)      = false
isvcrelax(calculation::DFCalculation)  = false
isrelax(calculation::DFCalculation)    = false
isprojwfc(calculation::DFCalculation)  = false
ismagnetic(calculation::DFCalculation) = false
issoc(calculation::DFCalculation)      = false

# QE interface
isbands(c::DFCalculation{QE})   = c[:calculation] == "bands"
isnscf(c::DFCalculation{QE})    = c[:calculation] == "nscf"
isscf(c::DFCalculation{QE})     = c[:calculation] == "scf"
isvcrelax(c::DFCalculation{QE}) = c[:calculation] == "vc-relax"
isrelax(c::DFCalculation{QE})   = c[:calculation] == "relax"

function ispw(c::DFCalculation{QE})
    return isbands(c) || isnscf(c) || isscf(c) || isvcrelax(c) || isrelax(c)
end

isprojwfc(c::DFCalculation{QE}) = DFC.hasexec(c, "projwfc.x")
ishp(c::DFCalculation{QE})      = DFC.hasexec(c, "hp.x")
issoc(c::DFCalculation{QE})     = flag(c, :lspinorb) == true

function ismagnetic(c::DFCalculation{QE})
    return (hasflag(c, :nspin) && c[:nspin] > 0.0) ||
           (hasflag(c, :total_magnetization) && c[:total_magnetization] != 0.0)
end


# ELK

infilename(calculation::DFCalculation{Elk}) = "elk.in"

isbandscalc(calculation::DFCalculation{Elk}) = calculation.name == "20"

isnscfcalc(calculation::DFCalculation{Elk}) = calculation.name == "elk2wannier" #nscf == elk2wan??

isscfcalc(calculation::DFCalculation{Elk}) = calculation.name ∈ ["0", "1"]


"""
    set_name!(calculation::DFCalculation, name::AbstractString)

Sets `calculation.name`, and `calculation.infile` and `calculation.outfile` to conform
with the new `name`.
"""
function set_name!(calculation::DFCalculation, name::AbstractString; print = true)
    calculation.name = name
    calculation.infile = name * splitext(infilename(calculation))[2]
    calculation.outfile = name * splitext(outfilename(calculation))[2]
    print &&
        @info "\ncalculation.name = $name\ncalculation.infile = $(infilename(calculation))\ncalculation.outfile = $(outfilename(calculation))"
    return name
end

#
# Directory interface
#
DFC.Utils.searchdir(i::DFCalculation, glob) = joinpath.((i,), searchdir(dir(i), glob))

"""
    joinpath(calc::DFCalculation, path...) = joinpath(dir(calc), path...)
"""
Base.joinpath(c::DFCalculation, path...) = joinpath(dir(c), path...)

#
# Flag Dict interace
#

"""
    getindex(calculation::DFCalculation, n::Symbol)

Returns the flag with given symbol.

    getindex(job::DFJob, flag::Symbol)
    
Searches through the job's calculations for the requested flag.
A `Dict` will be returned with calculation names as the keys and the flags as values.
"""
Base.getindex(c::DFCalculation, n::Symbol) =
    get(c, n, throw(KeyError(n)))

Base.haskey(c::DFCalculation, n::Symbol) = haskey(c.flags, n)
Base.get(c::DFCalculation, args...) = get(flags(c), args...)

hasflag(c::DFCalculation, s::Symbol) = haskey(flags(c), s)
hasflag(c::DFCalculation, s) = false

"""
    setindex!(calculation::DFCalculation, value, flag::Symbol)

Sets flags.

    setindex!(job::DFJob, value, flag::Symbol)

Set `flag` in all the appropriate calculations to the `value`.
"""
Base.setindex!(calculation::DFCalculation, dat, key) = set_flags!(calculation, key => dat)

function Base.setindex!(job::DFJob, value, key::Symbol)
    for calculation in calculations(job)
        calculation[key] = value
    end
end

"""
    set_flags!(calculation::DFCalculation, flags::Pair{Symbol, Any}...; print=true)

Sets multiple flags in one go. Flag validity and type are verified.

    set_flags!(job::DFJob, calculations::Vector{<:DFCalculation}, flags::Pair{Symbol,<:Any}...; print=true)
    set_flags!(job::DFJob, flags::Pair{Symbol,<:Any}...; print=true)

Sets the flags in the names to the flags specified.
This only happens if the specified flags are valid for the names.

The values that are supplied will be checked whether they are valid.
"""
function set_flags!(calculation::DFCalculation{T}, flags...; print = true) where {T}
    found_keys = Symbol[]
    for (flag, value) in flags
        flag_type = flagtype(calculation, flag)
        if flag_type !== Nothing
            !(flag in found_keys) && push!(found_keys, flag)
            try
                if isa(value, AbstractVector{<:AbstractVector}) &&
                   flag_type <: AbstractVector
                    value = [convert.(eltype(flag_type), v) for v in value]
                else
                    value = convert(flag_type, value)
                end
            catch
                print &&
                    (@warn "Filename '$(name(calculation))':\n  Could not convert '$value' into '$flag_type'.\n    Flag '$flag' not set.\n")
                continue
            end
            old_data = haskey(calculation.flags, flag) ? calculation.flags[flag] : ""
            calculation.flags[flag] = value
            print &&
                (@info "$(name(calculation)):\n  -> $flag:\n      $old_data set to: $value\n")
        else
            print &&
                @warn "Flag $flag was ignored since it could not be found in the allowed flags for calculation $(name(calculation))."
        end
    end
    return found_keys, calculation
end

function set_flags!(job::DFJob, flags...; kwargs...)
    return set_flags!(job, calculations(job), flags...; kwargs...)
end

"""
    rm_flags!(calculation::DFCalculation, flags::Symbol...)
    rm_flags!(job::DFJob, flags::Symbol...)

Remove the specified flags.
"""
function rm_flags!(calculation::DFCalculation, flags::Symbol...; print = true)
    for flag in flags
        if haskey(calculation.flags, flag)
            pop!(calculation.flags, flag, false)
            print && (@info "Removed flag '$flag' from calculation '$(name(calculation))'")
        end
    end
    return calculation
end

function rm_flags!(job::DFJob, flags::Symbol...; kwargs...)
    for c in calculations(job)
        rm_flags!(c, flags...; kwargs...)
    end
end



function Base.:(==)(d1::InputData, d2::InputData)
    return all(x -> getfield(d1, x) == getfield(d2, x), fieldnames(InputData))
end

package(::DFCalculation{T}) where {T} = T

data(calculation::DFCalculation) = calculation.data

hasexec(calculation::DFCalculation, ex::AbstractString) = exec(calculation, ex) != nothing

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

ψ_cutoff_flag(c::DFCalculation) = nothing
ρ_cutoff_flag(c::DFCalculation) = nothing

include("qe/calculation.jl")
include("wannier90/calculation.jl")
include("elk/calculation.jl")
