#these are all the control data, they hold the flags that guide the calculation

"""
    InputData(name::Symbol, option::Symbol, data::Any)

Represents a more structured block of input data.
e.g. `InputData(:k_points, :automatic, (6,6,6,1,1,1))`
would be translated for a QE calculation into
```
K_POINTS(automatic)
6 6 6 1 1 1
```
"""
mutable struct InputData
    name   ::Symbol
    option ::Symbol
    data   ::Any
end
name(data::InputData) = data.name
Base.:(==)(d1::InputData, d2::InputData) =
    all(x -> getfield(d1, x) == getfield(d2, x), fieldnames(InputData))

"""
    DFCalculation{P<:Package}(name::String;
                              dir      ::String = "",
                              flags    ::AbstractDict = Dict{Symbol, Any}(),
                              data     ::Vector{InputData} = InputData[],
                              execs    ::Vector{Exec},
                              run      ::Bool = true,
                              outdata  ::AbstractDict = Dict{Symbol, Any}(),
                              infile   ::String = P == Wannier90 ? name * ".win" : name * ".in",
                              outfile  ::String = name * ".out")

A full representation of a *DFT* calculation of package `P`,
holding the `flags` that will be written to the `infile`,
the executables in `execs` and the output written by the calculation to the `outfile`.
`outdata` stores the parsed calculation output after it was read at least once.
The `run` field indicates whether the calculation should be actually performed,
e.g. if `run=false` the corresponding line will be commented out in the job script.

    DFCalculation{P<:Package}(name::AbstractString, flags::Pair{Symbol, Any}...; kwargs...)

Create a `DFCalculation` from `name` and `flags`, other `kwargs...` will be passed to the constructor.

    DFCalculation(template::DFCalculation, name::AbstractString, flags::Pair{Symbol, Any}...; excs=deepcopy(execs(template)), run=true, data=nothing, dir=copy(template.dir))

Creates a new `DFCalculation` from the `template`, setting the `flags` of the newly created one to the specified ones.
"""
@with_kw_noshow mutable struct DFCalculation{P <: Package}
    name     ::String
    dir      ::String = ""
    flags    ::AbstractDict = SymAnyDict()
    data     ::Vector{InputData} = InputData[]
    execs    ::Vector{Exec}
    run      ::Bool = true
    outdata  ::SymAnyDict=SymAnyDict()
    infile   ::String = P == Wannier90 ? name * ".win" : name * ".in"
    outfile  ::String = name * ".out"
    function DFCalculation{P}(name, dir, flags, data, execs, run, outdata, infile, outfile) where {P <: Package}
        out = new{P}(name, dir, SymAnyDict(), data, execs, run, outdata, infile, outfile)
        set_flags!(out, flags...; print=false)
        for (f, v) in flags
            if !hasflag(out, f)
                @warn "Flag $f was not valid for calculation $name."
            end
        end
        return out
    end
end
DFCalculation{P}(name, dir, flags, data, execs, run) where P<:Package = DFCalculation{P}(name, abspath(dir), flags, data, execs, run, SymAnyDict(), P == Wannier90 ? name * ".win" : name * ".in", P == Wannier90 ? name * ".wout" : name * ".out")

DFCalculation{P}(name, flags...; kwargs...) where P<:Package =
    DFCalculation{P}(name=name, flags=flags; kwargs...)

function DFCalculation(template::DFCalculation, name, newflags...;
                 excs = deepcopy(execs(template)),
                 run  = true,
                 data = nothing,
                 dir  = deepcopy(template.dir))
    newflags = Dict(newflags...)

    calculation          = deepcopy(template)
    calculation.name     = name
    calculation.execs    = excs
    calculation.run      = run
    calculation.dir      = dir
    set_flags!(calculation, newflags..., print=false)

    if data != nothing
        for (name, (option, data)) in data
            set_data!(calculation, name, data, option=option, print=false)
        end
    end
    return calculation
end

name(calculation::DFCalculation)  = calculation.name
dir(calculation::DFCalculation)   = calculation.dir
flags(calculation::DFCalculation) = calculation.flags
set_dir!(calculation::DFCalculation, dir) = (calculation.dir = dir)
name_ext(calculation::DFCalculation, ext)   = name(calculation) * ext
infilename(calculation::DFCalculation)      = calculation.infile
outfilename(calculation::DFCalculation)     = calculation.outfile
inpath(calculation::DFCalculation)          = joinpath(calculation,  infilename(calculation))
outpath(calculation::DFCalculation)         = joinpath(calculation,  outfilename(calculation))


"""
    joinpath(calc::DFCalculation, path...)
"""
Base.joinpath(c::DFCalculation, path...) = joinpath(dir(c), path...)
hasflag(c::DFCalculation, s::Symbol) = haskey(flags(c), s)
hasflag(c::DFCalculation, s) = false

function flag(calculation::DFCalculation, flag::Symbol)
    if hasflag(calculation, flag)
        return calculation.flags[flag]
    end
end

Base.eltype(::DFCalculation{P}) where P = P
package(::DFCalculation{P}) where P = P

data(calculation::DFCalculation)  = calculation.data

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
    convert_flags!(calculation)
end

function setoradd!(datas::Vector{InputData}, data::InputData)
    found = false
    for (i, d) in enumerate(datas)
        if d.name == data.name
            datas[i] = data
            found = true
            break
        end
    end
    if !found
        push!(datas, data)
    end
end

"""
    set_data!(calculation::DFCalculation, data::InputData)

Adds the given data to the calculation. Should put it in the correct arrays.
"""
function set_data!(calculation::DFCalculation, data::InputData)
    setoradd!(calculation.data, data)
    return calculation
end

isbands(calculation::DFCalculation)        = false
isnscf(calculation::DFCalculation)         = false
isscf(calculation::DFCalculation)          = false
isvcrelax(calculation::DFCalculation)      = false
isprojwfc(calculation::DFCalculation)      = false
ismagnetic(calculation::DFCalculation)    = false
issoc(calculation::DFCalculation) = false


#TODO review this!
outdata(calculation::DFCalculation) = calculation.outdata
hasoutput(calculation::DFCalculation) = !isempty(outdata(calculation)) || hasoutfile(calculation)

hasoutfile(calculation::DFCalculation) = ispath(outpath(calculation))
hasnewout(calculation::DFCalculation, time) = mtime(outpath(calculation)) > time

set_cutoffs!(calculation::DFCalculation, args...) = @warn "Setting cutoffs is not implemented for package $(package(calculation))"

function hasoutput_assert(i::DFCalculation)
    @assert hasoutfile(i) "Please specify an calculation that has an outputfile."
end

function hasexec_assert(i::DFCalculation, exec::String)
    @assert hasexec(i, exec) "Please specify an calculation with $exec as it's executable."
end

Base.:(==)(i1::DFCalculation, i2::DFCalculation) =
    all(x -> x == :outdata ? true : getfield(i1, x) == getfield(i2, x), fieldnames(DFCalculation))

searchdir(i::DFCalculation, glob) = joinpath.((i,), searchdir(dir(i), glob))

ψ_cutoff_flag(::DFCalculation{P}) where {P} = ψ_cutoff_flag(P)
ρ_cutoff_flag(::DFCalculation{P}) where {P} = ρ_cutoff_flag(P)

pdos(calculation::DFCalculation, args...) =
    @error "pdos reading not implemented for package $(package(calculation))."

Emin_from_projwfc(calculation::DFCalculation, args...) =
    @error "Emin_from_projwfc is not implemented for package $(package(calculation))."

include("qe/calculation.jl")
include("elk/calculation.jl")
include("wannier90/calculation.jl")



