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
infilename(input::DFInput{Elk})        = "elk.in"
outfilename(input::DFInput{QE})        = name_ext(input, ".out")
outfilename(input::DFInput{Wannier90}) = name_ext(input, ".wout")
outfilename(input::DFInput{Elk})       = "elk.out"
inpath(input::DFInput)                 = joinpath(dir(input),  infilename(input))
outpath(input::DFInput)                = joinpath(dir(input),  outfilename(input))

hasflag(i::DFInput, s::Symbol) = haskey(flags(i), s)

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
    if isvcrelaxcalc(input)
	    #this is to make sure &ions and &cell are there in the input 
	    !hasflag(input, :ion_dynamics)  && setflags!(input, :ion_dynamics  => "bfgs", print=false)
	    !hasflag(input, :cell_dynamics) && setflags!(input, :cell_dynamics => "bfgs", print=false)
    end
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

isbandscalc(input::DFInput{QE})    = flag(input, :calculation) == "bands"
isbandscalc(input::DFInput{Elk})   = input.name == "20"
isbandscalc(input::DFInput)        = false

isnscfcalc(input::DFInput{QE})     = flag(input, :calculation) == "nscf"
isnscfcalc(input::DFInput{Elk})    = input.name == "elk2wannier" #nscf == elk2wan??
isnscfcalc(input::DFInput)         = false

isscfcalc(input::DFInput{QE})      = flag(input, :calculation) == "scf"
isscfcalc(input::DFInput{Elk})     = input.name ∈ ["0", "1"]
isscfcalc(input::DFInput)          = false

iscolincalc(input::DFInput{QE})    = flag(input, :nspin) == 2
iscolincalc(input::DFInput)        = false

isnoncolincalc(input::DFInput{QE}) = flag(input, :noncolin) == true
isnoncolincalc(input::DFInput)     = false

isvcrelaxcalc(input::DFInput{QE})  = flag(input, :calculation) == "vc-relax"
isvcrelaxcalc(input::DFInput)      = false

isprojwfccalc(input::DFInput{QE})  = hasexec(input, "projwfc.x")
isprojwfccalc(input::DFInput)      = false

isdftucalc(input::DFInput{QE})     = flag(input, :lda_plus_u) == true
isdftucalc(input::DFInput)         = false

ismagneticcalc(input::DFInput{QE}) = flag(input, :nspin) ∈ [2, 4] || (flag(input, :lda_plus_u) == true && flag(input, :noncolin) == true)
ismagneticcalc(input::DFInput)     = false

#TODO review this!
outdata(input::DFInput) = input.outdata
hasoutput(input::DFInput) = !isempty(outdata(input))

hasoutfile(input::DFInput) = ispath(outpath(input))
hasnewout(input::DFInput, time) = mtime(outpath(input)) > time

readoutput(input::DFInput{QE}) = qe_read_output(input)
readoutput(input::DFInput{Wannier90}) = wan_read_output(outpath(input))

pseudodir(input::DFInput{QE}) = flag(input, :pseudo_dir)

setcutoffs!(input::DFInput, args...) = @warn "Setting cutoffs is not implemented for package $(package(input))"
setcutoffs!(input::DFInput{QE}, ecutwfc, ecutrho) = setflags!(input, :ecutwfc => ecutwfc, :ecutrho=>ecutrho)

function Emin_from_projwfc(structure::AbstractStructure, projwfc::DFInput{QE}, threshold::Number)
    hasoutput_assert(projwfc)
    hasexec_assert(projwfc, "projwfc.x")
    if !haskey(outdata(projwfc), :states)
        states, bands = qe_read_projwfc(outpath(projwfc))
        outdata(projwfc)[:states] = states
        outdata(projwfc)[:bands]  = bands
    else
        states, bands = outdata(projwfc)[:states], outdata(projwfc)[:bands]
    end

    mask = zeros(length(states))
    for (atid, at) in enumerate(atoms(structure))
        projs = projections(at)

        if isempty(projs)
            continue
        end

        stateids = Int[]
        for proj in projs
            orb = orbital(proj)
            push!.((stateids,), findall(x -> x.atom_id == atid && x.l == orb.l, states))
        end
        mask[stateids] .= 1.0

    end
    Emin = -10000.0
    for b in bands
        ψ = mean(b.extra[:ψ])
        tot_relevant_occupation = dot(mask, ψ)

        if tot_relevant_occupation > threshold
            Emin = minimum(b.eigvals)
            break
        end

    end
    if Emin == -10000.0
        error("Couldn't find any band with occupation of relevant projections above $threshold.")
    end
    return Emin
end

#asserts
function iscalc_assert(i::DFInput{QE}, calc)
    @assert flag(i, :calculation) == calc "Please provide a valid '$calc' calculation."
end

function hasoutput_assert(i::DFInput)
    @assert hasoutfile(i) "Please specify an input that has an outputfile."
end

function hasexec_assert(i::DFInput, exec::String)
    @assert hasexec(i, exec) "Please specify an input with $exec as it's executable."
end


#TODO Temporary handling of HubbardU situation
function set_hubbard_flags!(input::DFInput{QE}, str::AbstractStructure{T}) where {T}
	if isdftucalc(input)
		u_ats = unique(atoms(str))
		setflags!(input,:Hubbard_U     => map(x -> dftu(x).U , u_ats); print=false)
		setflags!(input,:Hubbard_alpha => map(x -> dftu(x).α , u_ats); print=false)
		setflags!(input,:Hubbard_beta  => map(x -> dftu(x).β , u_ats); print=false)
		Jmap = map(x -> dftu(x).J, u_ats)
		setflags!(input, :Hubbard_J => reshape(collect(Iterators.flatten(Jmap)), length(u_ats), length(Jmap[1])); print=false)
		setflags!(input, :Hubbard_J0=> map(x -> dftu(x).J0, u_ats); print=false)
	end
end

function set_starting_magnetization_flags!(input::DFInput{QE}, str::AbstractStructure{T}) where {T}
	if ismagneticcalc(input)
		u_ats = unique(atoms(str))
		mags  = magnetization.(u_ats)
		starts= T[]
		θs    = T[]
		ϕs    = T[]
		for m in mags
			if norm(m) == 0
				push!.((starts, θs, ϕs), 0.0)
			else
				θ = acos(m[3])/norm(m) * 180/π
				ϕ = atan(m[2], m[1])   * 180/π
				start = 1
				push!(θs, θ)
				push!(ϕs, ϕ)
				push!(starts, start)
			end
		end
		setflags!(input, :starting_magnetization => starts; print=false)
		setflags!(input, :angle1 => θs; print=false)
		setflags!(input, :angle2 => ϕs; print=false)
	end
end








