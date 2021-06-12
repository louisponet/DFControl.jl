#these are all the control data, they hold the flags that guide the calculation
mutable struct InputData
    name   ::Symbol
    option ::Symbol
    data   ::Any
end
name(data::InputData) = data.name
Base.:(==)(d1::InputData, d2::InputData) =
    all(x -> getfield(d1, x) == getfield(d2, x), fieldnames(InputData))

@with_kw_noshow mutable struct DFInput{P <: Package}
    name     ::String
    dir      ::String = ""
    flags    ::AbstractDict = SymAnyDict()
    data     ::Vector{InputData} = InputData[]
    execs    ::Vector{Exec}
    run      ::Bool = true
    outdata  ::SymAnyDict=SymAnyDict()
    infile   ::String = P == Wannier90 ? name * ".win" : name * ".in"
    outfile  ::String = name * ".out"
end
DFInput{P}(name, dir, flags, data, execs, run) where P<:Package = DFInput{P}(name, abspath(dir), flags, data, execs, run, SymAnyDict(), P == Wannier90 ? name * ".win" : name * ".in", P == Wannier90 ? name * ".wout" : name * ".out")

"""
    DFInput(template::DFInput, name::AbstractString, flags::Pair{Symbol}...; excs=deepcopy(execs(template)), run=true, data=nothing, dir=copy(template.dir))

Creates a new `DFInput` from a template `DFInput`, setting the newflags in the new one.
"""
function DFInput(template::DFInput, name, newflags...; excs=deepcopy(execs(template)), run=true, data=nothing, dir=deepcopy(template.dir))
    newflags = Dict(newflags...)

    input          = deepcopy(template)
    input.name     = name
    input.execs    = excs
    input.run      = run
    input.dir      = dir
    set_flags!(input, newflags..., print=false)

    if data != nothing
        for (name, (option, data)) in data
            set_data!(input, name, data, option=option, print=false)
        end
    end
    return input
end

"""
    DFInput{P}(name::AbstractString, execs::Vector{Exec}, flags::Pair{Symbol}...; kwargs...) where {P<:Package}

Create a input from the name and flags, other fields will set to default values.
"""
function DFInput{P}(name::AbstractString, execs::Vector{Exec}, flags::Pair{Symbol}...; kwargs...) where {P<:Package}
    out = DFInput{P}(name=name, execs=execs; kwargs...)
    set_flags!(out, flags..., print=false)
    return out
end

name(input::DFInput)  = input.name
dir(input::DFInput)   = input.dir
Base.joinpath(i::DFInput, path) = joinpath(dir(i), path)
flags(input::DFInput) = input.flags
set_dir!(input::DFInput, dir) = (input.dir = dir)
name_ext(input::DFInput, ext)   = name(input) * ext
infilename(input::DFInput)      = input.infile
outfilename(input::DFInput)     = input.outfile
inpath(input::DFInput)          = joinpath(input,  infilename(input))
outpath(input::DFInput)         = joinpath(input,  outfilename(input))

hasflag(i::DFInput, s::Symbol) = haskey(flags(i), s)
hasflag(i::DFInput, s) = false

function flag(input::DFInput, flag::Symbol)
    if hasflag(input, flag)
        return input.flags[flag]
    end
end

Base.eltype(::DFInput{P}) where P = P
package(::DFInput{P}) where P = P

data(input::DFInput)  = input.data

execs(input::DFInput) = input.execs
hasexec(input::DFInput, ex::AbstractString) = exec(input, ex) != nothing
set_flow!(input::DFInput, run) = input.run = run

"Runs through all the set flags and checks if they are allowed and set to the correct value"
function convert_flags!(input::DFInput)
    for (flag, value) in flags(input)
        flagtype_ = flagtype(input, flag)
        if flagtype_ == Nothing
            @warn "Flag $flag was not found in allowed flags for exec $(execs(input)[2]). Removing flag."
            rm_flags!(input, flag)
            continue
        end
        if !(isa(value, flagtype_) || eltype(value) <: flagtype_)
            try
                if isbitstype(eltype(value))
                    if length(value) > 1
                        flags(input)[flag] = convert(flagtype_, value)
                    else
                        flags(input)[flag] = convert(eltype(flagtype_), value)
                    end
                else
                    flags(input)[flag] = convert.(flagtype_, value)
                end
            catch
                error("Input $(name(input)): Could not convert :$flag of value $value to the correct type ($flagtype_), please set it to the correct type.")
            end
        end
    end
end

#TODO implement abinit and wannier90
"""
    sanitize_flags!(input::DFInput, str::AbstractStructure)

Cleans up flags, i.e. remove flags that are not allowed and convert all
flags to the correct types.
Tries to correct common errors for different input types.
"""
function sanitize_flags!(input::DFInput, str::AbstractStructure)
    convert_flags!(input)
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
    set_data!(input::DFInput, data::InputData)

Adds the given data to the input. Should put it in the correct arrays.
"""
function set_data!(input::DFInput, data::InputData)
    setoradd!(input.data, data)
    return input
end

isbandscalc(input::DFInput)        = false
isnscfcalc(input::DFInput)         = false
isscfcalc(input::DFInput)          = false
isvcrelaxcalc(input::DFInput)      = false
isprojwfccalc(input::DFInput)      = false

issoccalc(input::DFInput) = false


#TODO review this!
outdata(input::DFInput) = input.outdata
hasoutput(input::DFInput) = !isempty(outdata(input)) || hasoutfile(input)

hasoutfile(input::DFInput) = ispath(outpath(input))
hasnewout(input::DFInput, time) = mtime(outpath(input)) > time

set_cutoffs!(input::DFInput, args...) = @warn "Setting cutoffs is not implemented for package $(package(input))"

function hasoutput_assert(i::DFInput)
    @assert hasoutfile(i) "Please specify an input that has an outputfile."
end

function hasexec_assert(i::DFInput, exec::String)
    @assert hasexec(i, exec) "Please specify an input with $exec as it's executable."
end

Base.:(==)(i1::DFInput, i2::DFInput) =
    all(x -> x == :outdata ? true : getfield(i1, x) == getfield(i2, x), fieldnames(DFInput))

searchdir(i::DFInput, glob) = joinpath.((i,), searchdir(dir(i), glob))

ψ_cutoff_flag(::DFInput{P}) where {P} = ψ_cutoff_flag(P)
ρ_cutoff_flag(::DFInput{P}) where {P} = ρ_cutoff_flag(P)


include("qe/input.jl")
include("elk/input.jl")
include("wannier90/input.jl")



