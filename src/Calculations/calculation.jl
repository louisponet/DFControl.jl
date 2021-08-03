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
    name   :: Symbol
    option :: Symbol
    data   :: Any
end

StructTypes.StructType(::Type{InputData}) = StructTypes.Struct()


abstract type Package end
struct NoPackage <: Package end
struct Wannier90 <: Package end
struct QE <: Package end
struct Abinit <: Package end
struct Elk <: Package end
StructTypes.StructType(::Type{<:Package}) = StructTypes.Struct()
export Wannier90, QE, Abinit, Elk


"""
    Calculation{P<:Package}(name    ::String;
                              dir     ::String = "",
                              flags   ::AbstractDict = Dict{Symbol, Any}(),
                              data    ::Vector{InputData} = InputData[],
                              execs   ::Vector{Exec},
                              run     ::Bool = true,
                              outdata ::AbstractDict = Dict{Symbol, Any}(),
                              infile  ::String = P == Wannier90 ? name * ".win" : name * ".in",
                              outfile ::String = name * ".out")

The representation of a *DFT* calculation of package `P`,
holding the `flags` that will be written to the `infile`,
the executables in `execs` and the output written by the calculation to the `outfile`.
It essentially represents a line in a job script similar to `exec1 exec2 < infile.in > outfile.out`. 
`outdata` stores the parsed calculation output after it was read at least once.
The `run` field indicates whether the calculation should be actually performed,
e.g. if `run=false` the corresponding line will be commented out in the job script.

    Calculation{P<:Package}(name::AbstractString, flags::Pair{Symbol, Any}...; kwargs...)

Create a [`Calculation`](@ref) from `name` and `flags`, other `kwargs...` will be passed to the constructor.

    Calculation(template::Calculation, name::AbstractString, flags::Pair{Symbol, Any}...; excs=deepcopy(execs(template)), run=true, data=nothing, dir=copy(template.dir))

Creates a new [`Calculation`](@ref) from the `template`, setting the `flags` of the newly created one to the specified ones.
"""
@with_kw_noshow mutable struct Calculation{P<:Package}
    name::String
    dir::String = ""
    flags::SymAnyDict = SymAnyDict()
    data::Vector{InputData} = InputData[]
    execs::Vector{Exec}
    run::Bool = true
    outdata::SymAnyDict = SymAnyDict()
    infile::String = P == Wannier90 ? name * ".win" : name * ".in"
    outfile::String = P == Wannier90 ? name * ".wout" : name * ".out"
    function Calculation{P}(name, dir, flags, data, execs, run, outdata, infile,
                              outfile) where {P<:Package}
        out = new{P}(name, dir, SymAnyDict(), data, execs, run, outdata, infile, outfile)
        set_flags!(out, flags...; print = false)
        for (f, v) in flags
            if !hasflag(out, f)
                @warn "Flag $f was not valid for calculation $name."
            end
        end
        return out
    end
end
function Calculation{P}(name, dir, flags, data, execs, run) where {P<:Package}
    return Calculation{P}(name, abspath(dir), flags, data, execs, run, SymAnyDict(),
                            P == Wannier90 ? name * ".win" : name * ".in",
                            P == Wannier90 ? name * ".wout" : name * ".out")
end

function Calculation{P}(name, flags...; kwargs...) where {P<:Package}
    return Calculation{P}(; name = name, flags = flags, kwargs...)
end

function Calculation(template::Calculation, name, newflags...;
                       excs = deepcopy(execs(template)), run  = true, data = nothing,
                       dir  = deepcopy(template.dir))
    newflags = Dict(newflags...)

    calculation       = deepcopy(template)
    calculation.name  = name
    calculation.execs = excs
    calculation.run   = run
    calculation.dir   = dir
    set_flags!(calculation, newflags...; print = false)

    if data != nothing
        for (name, (option, data)) in data
            set_data!(calculation, name, data; option = option, print = false)
        end
    end
    return calculation
end

# Calculation() = Calculation{NoPackage}(package=NoPackage())
StructTypes.StructType(::Type{<:Calculation}) = StructTypes.Mutable()

set_dir!(c::Calculation, dir) = (c.dir = dir)
inpath(c::Calculation)        = joinpath(c, c.infile)
outpath(c::Calculation)       = joinpath(c, c.outfile)

# Interface Functions
isbands(c::Calculation)    = false
isnscf(c::Calculation)     = false
isscf(c::Calculation)      = false
isvcrelax(c::Calculation)  = false
isrelax(c::Calculation)    = false
isprojwfc(c::Calculation)  = false
ismagnetic(c::Calculation) = false
issoc(c::Calculation)      = false
outfiles(c::Calculation)   = [outpath(c)]

# QE interface
isbands(c::Calculation{QE})   = c[:calculation] == "bands"
isnscf(c::Calculation{QE})    = c[:calculation] == "nscf"
isscf(c::Calculation{QE})     = c[:calculation] == "scf"
isvcrelax(c::Calculation{QE}) = c[:calculation] == "vc-relax"
isrelax(c::Calculation{QE})   = c[:calculation] == "relax"

function ispw(c::Calculation{QE})
    return isbands(c) || isnscf(c) || isscf(c) || isvcrelax(c) || isrelax(c)
end

issoc(c::Calculation{QE})     = c[:lspinorb] == true

function ismagnetic(c::Calculation{QE})
    return (hasflag(c, :nspin) && c[:nspin] > 0.0) ||
           (hasflag(c, :total_magnetization) && c[:total_magnetization] != 0.0)
end

function outfiles(c::Calculation{QE})
    files = [outpath(c)]
    for (is, fuzzies) in zip(("projwfc.x", "hp.x"), (("pdos",), ("Hubbard_parameters",)))
        if any(x -> x.exec == is, c.execs)
            for f in fuzzies
                append!(files, searchdir(c, f))
            end
        end
    end
    return filter(ispath, unique(files))
end

# ELK
infilename(c::Calculation{Elk}) = "elk.in"
isbandscalc(c::Calculation{Elk}) = c.name == "20"
isnscfcalc(c::Calculation{Elk}) = c.name == "elk2wannier" #nscf == elk2wan??
isscfcalc(c::Calculation{Elk}) = c.name ∈ ["0", "1"]

"""
    set_name!(c::Calculation, name::AbstractString)

Sets `calculation.name`, and `calculation.infile` and `calculation.outfile` to conform
with the new `name`.
"""
function set_name!(c::Calculation, name::AbstractString; print = true)
    c.name = name
    c.infile = name * splitext(c.infile)[2]
    c.outfile = name * splitext(c.outfile)[2]
    print &&
        @info "\nname = $name\ninfile = $(c.infile))\noutfile = $(c.outfile)"
    return name
end

#
# Directory interface
#
DFC.Utils.searchdir(i::Calculation, glob) = joinpath.((i,), searchdir(i.dir, glob))

"""
    joinpath(calc::Calculation, path...) = joinpath(calc.dir, path...)
"""
Base.joinpath(c::Calculation, path...) = joinpath(c.dir, path...)

Base.eltype(::Calculation{T}) where {T} = T

#
# Flag Dict interace
#

"""
    getindex(c::Calculation, n::Symbol)

Returns the flag with given symbol.

    getindex(job::DFJob, flag::Symbol)
    
Searches through the job's calculations for the requested flag.
A `Dict` will be returned with calculation names as the keys and the flags as values.
"""
Base.getindex(c::Calculation, n::Symbol) =
    get(c, n, throw(KeyError(n)))

Base.haskey(c::Calculation, n::Symbol) = haskey(c.flags, n)
Base.get(c::Calculation, args...) = get(c.flags, args...)
Base.pop!(c::Calculation, args...) = pop!(c.flags, args...)

hasflag(c::Calculation, s::Symbol) = haskey(c.flags, s)
hasflag(c::Calculation, s) = false

"""
    setindex!(c::Calculation, value, flag::Symbol)

Sets flags.

    setindex!(job::DFJob, value, flag::Symbol)

Set `flag` in all the appropriate calculations to the `value`.
"""
Base.setindex!(c::Calculation, dat, key) = set_flags!(c, key => dat)

"""
    set_flags!(c::Calculation, flags::Pair{Symbol, Any}...; print=true)

Sets multiple flags in one go. Flag validity and type are verified.

    set_flags!(job::DFJob, calculations::Vector{<:Calculation}, flags::Pair{Symbol,<:Any}...; print=true)
    set_flags!(job::DFJob, flags::Pair{Symbol,<:Any}...; print=true)

Sets the flags in the names to the flags specified.
This only happens if the specified flags are valid for the names.

The values that are supplied will be checked whether they are valid.
"""
function set_flags!(c::Calculation{T}, flags...; print = true) where {T}
    found_keys = Symbol[]
    for (flag, value) in flags
        flag_type = flagtype(c, flag)
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
                    (@warn "Filename '$(c.name)':\n  Could not convert '$value' into '$flag_type'.\n    Flag '$flag' not set.\n")
                continue
            end
            old_data = haskey(c.flags, flag) ? c.flags[flag] : ""
            c.flags[flag] = value
            print &&
                (@info "$(c.name):\  -> $flag:\n      $old_data set to: $value\n")
        else
            print &&
                @warn "Flag $flag was ignored since it could not be found in the allowed flags for calculation $(c.name)"
        end
    end
    return found_keys, c
end

"Runs through all the set flags and checks if they are allowed and set to the correct value"
function convert_flags!(calculation::Calculation)
    for (flag, value) in calculation.flags
        flagtype_ = DFC.flagtype(calculation, flag)
        if flagtype_ == Nothing
            @warn "Flag $flag was not found in allowed flags for exec $(execs(calculation)[2]). Removing flag."
            rm_flags!(calculation, flag)
            continue
        end
        if !(isa(value, flagtype_) || eltype(value) <: flagtype_)
            try
                if isbitstype(eltype(value))
                    if length(value) > 1
                        calculation[flag] = convert(flagtype_, value)
                    else
                        calculation[flag] = convert(eltype(flagtype_), value)
                    end
                else
                    calculation[flag] = convert.(flagtype_, value)
                end
            catch
                error("Input $(calculation.name): Could not convert :$flag of value $value to the correct type ($flagtype_), please set it to the correct type.")
            end
        end
    end
end

function Base.:(==)(d1::InputData, d2::InputData)
    return all(x -> getfield(d1, x) == getfield(d2, x), fieldnames(InputData))
end

function Base.:(==)(i1::Calculation, i2::Calculation)
    return all(x -> x in (:outdata, :dir, :run) ? true : getfield(i1, x) == getfield(i2, x),
               fieldnames(Calculation))
end

ψ_cutoff_flag(::Calculation) = nothing
ρ_cutoff_flag(::Calculation) = nothing
ψ_cutoff_flag(::Calculation{QE}) = :ecutwfc
ρ_cutoff_flag(::Calculation{Wannier90}) = :ecutrho

function set_cutoffs!(c::Calculation, ecutwfc, ecutrho)
    return set_flags!(c, ψ_cutoff_flag(c) => ecutwfc, ρ_cutoff_flag(c) => ecutrho)
end

function iscalc_assert(i::Calculation{QE}, calc)
    @assert flag(i, :calculation) == calc "Please provide a valid '$calc' calculation."
end

for f in (:cp, :mv)
    @eval function Base.$f(i::Calculation{QE}, dest::String; kwargs...)
        $f(inpath(i), joinpath(dest, i.infile); kwargs...)
        if hasoutfile(i)
            $f(outpath(i), joinpath(dest, i.outfile); kwargs...)
        end
        if any(x -> x.exec == "projwfc.x", c.execs)
            for f in searchdir(i, "pdos")
                $f(f, joinpath(dest, splitdir(f)[end]); kwargs...)
            end
        elseif any(x -> x.exec == "hp.x", c.execs)
            for f in searchdir(i, "Hubbard_parameters")
                $f(f, joinpath(dest, splitdir(f)[end]); kwargs...)
            end
        end
        #TODO add ph.x outfiles
    end
end

