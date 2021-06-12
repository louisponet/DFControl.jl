set_noncolin_flags!(i::DFInput{QE}) = set_flags!(i, :npsin => 4, :noncolin => true; print=false)

function sanitize_flags!(input::DFInput{QE}, structure::AbstractStructure)
    if isvcrelaxcalc(input)
	    #this is to make sure &ions and &cell are there in the input 
	    !hasflag(input, :ion_dynamics)  && set_flags!(input, :ion_dynamics  => "bfgs", print=false)
	    !hasflag(input, :cell_dynamics) && set_flags!(input, :cell_dynamics => "bfgs", print=false)
    end
    #TODO add all the required flags
    if exec(input, "pw.x") !== nothing
        @assert hasflag(input, :calculation) "Please set the flag for calculation with name: $(name(input))"
    end
    # setting hubbard and magnetization flags
    set_hubbard_flags!(input, structure)
    set_starting_magnetization_flags!(input, structure)
    
    # setting hubbard flags 
    pseudo_dir = pseudo(atoms(structure)[1]).dir # Pseudos should be all sanitized by now
       set_flags!(input, :pseudo_dir => pseudo_dir; print=false)
       
    convert_flags!(input)
end

isbandscalc(input::DFInput{QE})    = flag(input, :calculation) == "bands"
isnscfcalc(input::DFInput{QE})     = flag(input, :calculation) == "nscf"
isscfcalc(input::DFInput{QE})      = flag(input, :calculation) == "scf"
isvcrelaxcalc(input::DFInput{QE})  = flag(input, :calculation) == "vc-relax"
isprojwfccalc(input::DFInput{QE})  = hasexec(input, "projwfc.x")
ishpcalc(input::DFInput{QE}) = hasexec(input, "hp.x")
issoccalc(input::DFInput{QE}) = flag(input, :lspinorb) == true

readoutput(input::DFInput{QE}) = qe_read_output(input)

pseudodir(input::DFInput{QE}) = flag(input, :pseudo_dir)

set_cutoffs!(input::DFInput{QE}, ecutwfc, ecutrho) = set_flags!(input, :ecutwfc => ecutwfc, :ecutrho=>ecutrho)

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

function set_hubbard_flags!(input::DFInput{QE}, str::AbstractStructure{T}) where {T}
	u_ats = unique(atoms(str))
    isdftucalc = any(x -> dftu(x).U != 0 || dftu(x).J0 != 0.0 || sum(dftu(x).J) != 0 || sum(dftu(x).α) != 0, u_ats) || hasflag(input, :Hubbard_parameters)
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
    	set_flags!(input,
    	          :lda_plus_u    => true,
    	          :Hubbard_U     => map(x -> dftu(x).U, u_ats),
    	          :Hubbard_alpha => map(x -> dftu(x).α , u_ats),
    	          :Hubbard_beta  => map(x -> dftu(x).β , u_ats),
    	          :Hubbard_J     => Jarr,
    	          :Hubbard_J0    => map(x -> dftu(x).J0, u_ats);
	              print=false)
        isnc && set_flags!(input, :lda_plus_u_kind => 1; print=false)
    else
        rm_flags!(input, :lda_plus_u, :lda_plus_u_kind, :Hubbard_U, :Hubbard_alpha, :Hubbard_beta, :Hubbard_J, :Hubbard_J0, :U_projection_typel; print=false)
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
		set_flags!(input, :noncolin => true; print=false)
		rm_flags!(input, :nspin; print=false)
	elseif ismagcalc 
		for m in mags
			push!.((θs, ϕs), 0.0)
			if norm(m) == 0
				push!(starts, 0)
			else
				push!(starts, sign(sum(m))*norm(m))
			end
		end
		set_flags!(input, :nspin => 2; print=false)
	end
	set_flags!(input, :starting_magnetization => starts, :angle1 => θs, :angle2 => ϕs; print=false)
end

for f in (:cp, :mv)
    @eval function Base.$f(i::DFInput{QE}, dest::String; kwargs...)
        $f(inpath(i), joinpath(dest, infilename(i)); kwargs...)
        if hasoutfile(i)
            $f(outpath(i), joinpath(dest, outfilename(i)); kwargs...)
        end
        if isprojwfccalc(i)
            for f in searchdir(i, "pdos")
                $f(f, joinpath(dest, splitdir(f)[end]); kwargs...)
            end
        elseif ishpcalc(i)
            for f in searchdir(i, "Hubbard_parameters")
                $f(f, joinpath(dest, splitdir(f)[end]); kwargs...)
            end
        end
        #TODO add ph.x outfiles
    end
end
