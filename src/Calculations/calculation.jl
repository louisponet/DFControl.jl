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
    flags::Dict{Symbol,Any} = Dict{Symbol,Any}()
    data::Vector{InputData} = InputData[]
    execs::Vector{Exec}
    run::Bool = true
    outdata::Dict{Symbol,Any} = Dict{Symbol,Any}()
    infile::String = P == Wannier90 ? name * ".win" : name * ".in"
    outfile::String = P == Wannier90 ? name * ".wout" : name * ".out"
    function Calculation{P}(name, dir, flags, data, execs, run, outdata, infile,
                              outfile) where {P<:Package}
        out = new{P}(name, dir, Dict{Symbol,Any}(), data, execs, run, outdata, infile, outfile)
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
    return Calculation{P}(name, abspath(dir), flags, data, execs, run, Dict{Symbol,Any}(),
                            P == Wannier90 ? name * ".win" : name * ".in",
                            P == Wannier90 ? name * ".wout" : name * ".out")
end

function Calculation{P}(name, flags...; kwargs...) where {P<:Package}
    return Calculation{P}(; name = name, flags = flags, kwargs...)
end

function Calculation(template::Calculation, name, newflags...;
                       excs = deepcopy(template.execs), run  = true, data = nothing,
                       dir  = copy(template.dir))
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
ismagnetic(c::Calculation) = false
issoc(c::Calculation)      = false
outfiles(c::Calculation)   = [outpath(c)]


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

"""
    data(calculation::Calculation, n::Symbol)

The former returns `calculation.data`, the later -- the `InputData` with name `n`.
"""
data(calculation::Calculation, n::Symbol) = getfirst(x -> x.name == n, calculation.data)


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

function set_cutoffs!(c::Calculation, ecutwfc, ecutrho)
    return set_flags!(c, ψ_cutoff_flag(c) => ecutwfc, ρ_cutoff_flag(c) => ecutrho)
end


"""
    set_kpoints!(calculation::Calculation{QE}, k_grid::NTuple{3, Int}; print=true)
    set_kpoints!(calculation::Calculation{QE}, k_grid::NTuple{6, Int}; print=true)
    set_kpoints!(calculation::Calculation{QE}, k_grid::Vector{<:NTuple{4}}; print=true, k_option=:crystal_b)

Convenience function to set the `:k_points` data block of `calculation`.
The three different methods are targeted at `nscf`, `scf` or `vcrelax`,
and `bands` calculations, respectively.
For the `nscf` version an explicit list of `k_points` will be generated.

    set_kpoints!(calculation::Calculation{Wannier90}, k_grid::NTuple{3, Int})

Similar to the `nscf` targeted function in the sense that it will generate
an explicit list of `k_points`, adhering to the same rules as for the `nscf`.
The `mp_grid` flag will also automatically be set.
"""
function set_kpoints!(::Calculation{P}, args...; kwargs...) where {P}
    @error "set_kpoints! not implemented for package $P."
end

#-------- Generating new Calculations ---------- #

function calculation_from_kpoints(template::Calculation, newname, kpoints, newflags...)
    newcalc = Calculation(deepcopy(template); name = newname)
    set_flags!(newcalc, newflags...; print=false)
    set_name!(newcalc, newname)
    set_kpoints!(newcalc, kpoints; print = false)
    return newcalc
end

function gencalc_nscf(::Calculation{P}, args...) where {P}
    @error "gencalc_nscf is not implemented for package $P."
end

function gencalc_scf(::Calculation{P}, args...) where {P}
    @error "gencalc_scf is not implemented for package $P."
end

function gencalc_projwfc(::Calculation{P}, args...) where {P}
    @error "gencalc_projwfc is not implemented for package $P."
end

function gencalc_wan(::Calculation{P}, args...) where {P}
    @error "gencalc_wan is not implemented for package $P."
end

function gencalc_vcrelax(::Calculation{P}, args...) where {P}
    @error "gencalc_vcrelax is not implemented for package $P."
end

function gencalc_bands(::Calculation{P}, args...) where {P}
    @error "gencalc_bands is not implemented for package $P."
end

function isconverged(::Calculation{P}) where {P}
    @error "isconverged is not implemented for package $P."
end

"""
    kgrid(na, nb, nc, calculation)

Returns an array of k-grid points that are equally spaced, calculation can be either `:wan` or `:nscf`, the returned grids are appropriate as calculations for wannier90 or an nscf calculation respectively.
"""
kgrid(na, nb, nc, ::Calculation{T}) where {T} = kgrid(na, nb, nc, T)

