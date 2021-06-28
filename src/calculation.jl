#these are all the control data, they hold the flags that guide the calculation
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

outfiles(c::DFCalculation) = filter(ispath, [outpath(c)])

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


Base.eltype(::DFCalculation{P}) where {P} = P
package(::DFCalculation{P}) where {P} = P

data(calculation::DFCalculation) = calculation.data

execs(calculation::DFCalculation) = calculation.execs
hasexec(calculation::DFCalculation, ex::AbstractString) = exec(calculation, ex) != nothing
set_flow!(calculation::DFCalculation, run) = calculation.run = run

"Runs through all the set flags and checks if they are allowed and set to the correct value"
function convert_flags!(calculation::DFCalculation)
    for (flag, value) in flags(calculation)
        flagtype_ = flagtype(calculation, flag)
        if flagtype_ == Nothing
            @warn "Flag $flag was not found in allowed flags for exec $(execs(calculation)[2]). Removing flag."
            rm_flags!(calculation, flag)
            continue
        end
        if !(isa(value, flagtype_) || eltype(value) <: flagtype_)
            try
                if isbitstype(eltype(value))
                    if length(value) > 1
                        flags(calculation)[flag] = convert(flagtype_, value)
                    else
                        flags(calculation)[flag] = convert(eltype(flagtype_), value)
                    end
                else
                    flags(calculation)[flag] = convert.(flagtype_, value)
                end
            catch
                error("Input $(name(calculation)): Could not convert :$flag of value $value to the correct type ($flagtype_), please set it to the correct type.")
            end
        end
    end
end

#TODO implement abinit and wannier90
"""
    sanitize_flags!(calculation::DFCalculation, str::AbstractStructure)

Cleans up flags, i.e. remove flags that are not allowed and convert all
flags to the correct types.
Tries to correct common errors for different calculation types.
"""
function sanitize_flags!(calculation::DFCalculation, str::AbstractStructure)
    return convert_flags!(calculation)
end

"""
    set_data!(calculation::DFCalculation, data::InputData)

If an `InputData` with the same name as `data` is already in `calculation`, it will be overwritten. Otherwise `data` gets pushed to the list of `InputData` blocks.
"""
function set_data!(calculation::DFCalculation, data::InputData)
    id = findfirst(x -> x.name == data.name, calculation.data)
    if id === nothing
        push!(calculation.data, data)
    else
        calculation.data[id] = data
    end
    return calculation
end

isbands(calculation::DFCalculation)    = false
isnscf(calculation::DFCalculation)     = false
isscf(calculation::DFCalculation)      = false
isvcrelax(calculation::DFCalculation)  = false
isprojwfc(calculation::DFCalculation)  = false
ismagnetic(calculation::DFCalculation) = false
issoc(calculation::DFCalculation)      = false

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
    return all(x -> x == :outdata ? true : getfield(i1, x) == getfield(i2, x),
               fieldnames(DFCalculation))
end

searchdir(i::DFCalculation, glob) = joinpath.((i,), searchdir(dir(i), glob))

ψ_cutoff_flag(::DFCalculation{P}) where {P} = ψ_cutoff_flag(P)
ρ_cutoff_flag(::DFCalculation{P}) where {P} = ρ_cutoff_flag(P)

function pdos(calculation::DFCalculation, args...)
    @error "pdos reading not implemented for package $(package(calculation))."
end

function Emin_from_projwfc(calculation::DFCalculation, args...)
    @error "Emin_from_projwfc is not implemented for package $(package(calculation))."
end

function readoutput(c::DFCalculation, args...; kwargs...)
    @error "Output parsing for package $(package(c)) not implemented."
end

rm_outfiles(calc::DFCalculation) = rm.(outfiles(calc))


include("qe/calculation.jl")
include("elk/calculation.jl")
include("wannier90/calculation.jl")
