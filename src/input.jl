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
    DFInput(template::DFInput, name, newflags...; excs=deepcopy(execs(template)), run=true, data=nothing, dir=copy(template.dir))

Creates a new `DFInput` from a template `DFInput`, setting the newflags in the new one.
"""
function DFInput(template::DFInput, name, newflags...; excs=deepcopy(execs(template)), run=true, data=nothing, dir=deepcopy(template.dir))
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
infilename(input::DFInput)         = input.infile
infilename(input::DFInput{Elk})        = "elk.in"
outfilename(input::DFInput)        = input.outfile
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
    cleanflags!(input)
    if isvcrelaxcalc(input)
	    #this is to make sure &ions and &cell are there in the input 
	    !hasflag(input, :ion_dynamics)  && setflags!(input, :ion_dynamics  => "bfgs", print=false)
	    !hasflag(input, :cell_dynamics) && setflags!(input, :cell_dynamics => "bfgs", print=false)
    end
    #TODO add all the required flags
    if exec(input, "pw.x") !== nothing
        @assert hasflag(input, :calculation) "Please set the flag for calculation with name: $(name(input))"
    end
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

isvcrelaxcalc(input::DFInput{QE})  = flag(input, :calculation) == "vc-relax"
isvcrelaxcalc(input::DFInput)      = false

isprojwfccalc(input::DFInput{QE})  = hasexec(input, "projwfc.x")
isprojwfccalc(input::DFInput)      = false

issoccalc(input::DFInput{QE}) = flag(input, :lspinorb) == true
issoccalc(input::DFInput{Wannier90}) = flag(input, :spinors) == true
issoccalc(input::DFInput) = false

#TODO review this!
outdata(input::DFInput) = input.outdata
hasoutput(input::DFInput) = !isempty(outdata(input)) || hasoutfile(input)

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
        error("Couldn't find any band with occupation of relevant projections above $threshold, were any set in the structure?")
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
	u_ats = unique(atoms(str))
    isdftucalc = any(x -> dftu(x).U != 0 || dftu(x).J0 != 0.0 || sum(dftu(x).J) != 0, u_ats)
    isnc = isnoncolin(str)
	if isdftucalc
		Jmap = map(x -> copy(dftu(x).J), u_ats)
		Jdim = maximum(length.(Jmap))
		Jarr = zeros(Jdim, length(u_ats))
		for (i, J) in enumerate(Jmap)
			diff = Jdim - length(J)
			if diff > 0
				for d in 1:diff
					push!(J, zero(eltype(J)))
				end
			end
			Jarr[:, i] .= J
		end
    	setflags!(input,
    	          :lda_plus_u    => true,
    	          :Hubbard_U     => map(x -> dftu(x).U, u_ats),
    	          :Hubbard_alpha => map(x -> dftu(x).α , u_ats),
    	          :Hubbard_beta  => map(x -> dftu(x).β , u_ats),
    	          :Hubbard_J     => Jarr,
    	          :Hubbard_J0    => map(x -> dftu(x).J0, u_ats);
	              print=false)
        isnc && setflags!(input, :lda_plus_u_kind => 1; print=false)
	end
end

function set_starting_magnetization_flags!(input::DFInput{QE}, str::AbstractStructure{T}) where {T}
	u_ats = unique(atoms(str))
	mags  = magnetization.(u_ats)
	starts= T[]
	θs    = T[]
	ϕs    = T[]
    ismagcalc = ismagnetic(str)
    isnc =  isnoncolin(str)
	if (ismagcalc && isnc) || (flag(input, :noncolin) !== nothing && flag(input, :noncolin))
		for m in mags
			tm = normalize(m)
			if norm(m) == 0
				push!.((starts, θs, ϕs), 0.0)
			else
				θ = acos(tm[3]) * 180/π
				ϕ = atan(tm[2], tm[1])   * 180/π
				start = norm(m)
				push!(θs, θ)
				push!(ϕs, ϕ)
				push!(starts, start)
			end
		end
		setflags!(input, :noncolin => true; print=false)
		rmflags!(input, :nspin; print=false)
	elseif ismagcalc 
		for m in mags
			push!.((θs, ϕs), 0.0)
			if norm(m) == 0
				push!(starts, 0)
			else
				push!(starts, sign(sum(m))*norm(m))
			end
		end
		setflags!(input, :nspin => 2; print=false)
	end
	setflags!(input, :starting_magnetization => starts, :angle1 => θs, :angle2 => ϕs; print=false)
end

Base.:(==)(i1::DFInput, i2::DFInput) =
    all(x -> x == :outdata ? true : getfield(i1, x) == getfield(i2, x), fieldnames(DFInput)) 
