#these are all the control data, they hold the flags that guide the calculation
mutable struct InputData
    name   ::Symbol
    option ::Symbol
    data   ::Any
end
name(data::InputData) = data.name

@with_kw_noshow mutable struct DFInput{P <: Package}
    name     ::String
    dir      ::String
    flags    ::SymAnyDict
    data     ::Vector{InputData} = InputData[]
    execs    ::Vector{Exec}
    run      ::Bool
    outdata  ::SymAnyDict=SymAnyDict()
end
DFInput{P}(name, dir, flags, data, execs, run) where P<:Package = DFInput{P}(name, abspath(dir), flags, data, execs, run, SymAnyDict())

"""
    DFInput(template::DFInput, name, newflags...; runcommand=template.runcommand, run=true)

Creates a new `DFInput` from a template `DFInput`, setting the newflags in the new one.
"""
function DFInput(template::DFInput, name, newflags...; excs=execs(template), run=true, data=nothing, dir=template.dir)
    newflags = Dict(newflags...)

    input          = deepcopy(template)
    input.name     = name
    input.execs    = excs
    input.run      = run
    input.dir      = dir
    setflags!(input, newflags..., print=false)

    if data != nothing
        for (name, (option, data)) in data
            setdata!(input, name, data, option=option, print=false)
        end
    end
    return input
end

name(input::DFInput)  = input.name
dir(input::DFInput)   = input.dir
flags(input::DFInput) = input.flags
setdir!(input::DFInput, dir) = (input.dir = dir)
name_ext(input::DFInput, ext)          = name(input) * ext
infilename(input::DFInput{QE})         = name_ext(input, ".in")
infilename(input::DFInput{Wannier90})  = name_ext(input, ".win")
outfilename(input::DFInput{QE})        = name_ext(input, ".out")
outfilename(input::DFInput{Wannier90}) = name_ext(input, ".wout")
inpath(input::DFInput)                 = joinpath(dir(input),  infilename(input))
outpath(input::DFInput)                = joinpath(dir(input),  outfilename(input))

function flag(input::DFInput, flag::Symbol)
    if haskey(input.flags, flag)
        return input.flags[flag]
    end
end

Base.eltype(::DFInput{P}) where P = P
package(::DFInput{P}) where P = P

data(input::DFInput)  = input.data

execs(input::DFInput) = input.execs

setflow!(input::DFInput, run) = input.run = run

"Runs through all the set flags and checks if they are allowed and set to the correct value"
function cleanflags!(input::DFInput)
    for (flag, value) in flags(input)
        flagtype_ = flagtype(input, flag)
        if flagtype_ == Nothing
            @warn "Flag $flag was not found in allowed flags for exec $(execs(input)[2]). Removing flag."
            rmflags!(input, flag)
            continue
        end
        if !(eltype(value) <: flagtype_)
            try
                flags(input)[flag] = convert(flagtype_, value)
            catch
                error("Input $(name(input)): Could not convert :$flag of value $value to the correct type ($flagtype_), please set it to the correct type.")
            end
        end
    end
end

#TODO implement abinit and wannier90
"""
    sanitizeflags!(input::DFInput)

Cleans up flags, i.e. remove flags that are not allowed and convert all
flags to the correct types.
Tries to correct common errors for different input types.
"""
function sanitizeflags!(input::DFInput)
    cleanflags!(input)
end
function sanitizeflags!(input::DFInput{QE})
    setflags!(input, :outdir => "$(joinpath(dir(input), "outputs"))", print=false)
    cleanflags!(input)
    #TODO add all the required flags
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
    setdata!(input::DFInput, data::InputData)

Adds the given data to the input. Should put it in the correct arrays.
"""
function setdata!(input::DFInput, data::InputData)
    setoradd!(input.data, data)
    return input
end

isbandscalc(input::DFInput{QE}) = flag(input, :calculation) == "bands"
isnscfcalc(input::DFInput{QE}) = flag(input, :calculation) == "nscf"
isscfcalc(input::DFInput{QE}) = flag(input, :calculation) == "scf"
isspincalc(input::DFInput{QE}) = all(flag(input, :nspin) .!= [nothing, 1])

#TODO review this!
outdata(input::DFInput) = input.outdata
hasoutput(input::DFInput) = !isempty(outdata(input))

hasoutfile(input::DFInput) = ispath(outpath(input))
hasnewout(input::DFInput, time) = mtime(outpath(input)) > time

readoutput(input::DFInput{QE}) = qe_read_output(outpath(input))
readoutput(input::DFInput{Wannier90}) = read_wannier_output(outpath(input))

pseudodir(input::DFInput{QE}) = flag(input, :pseudo_dir)

setcutoffs!(input::DFInput, args...) = @warn "Setting cutoffs is not implemented for package $(package(input))"
setcutoffs!(input::DFInput{QE}, ecutwfc, ecutrho) = setflags!(input, :ecutwfc => ecutwfc, :ecutrho=>ecutrho)


function iscalc_assert(i::DFInput{QE}, calc)
    @assert flag(i, :calculation) == calc "Please provide a valid '$calc' calculation."
end

function hasoutput_assert(i::DFInput)
    @assert hasoutfile(i) "Please specify an input that has an outputfile."
end
