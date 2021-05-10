
#----------- Basic Interaction --------------#

"Returns the flag with given symbol."
function Base.getindex(input::DFInput, n::Symbol)
    if haskey(input.flags, n)
        return input.flags[n]
    else
        error("Flag:$n Not found in flags of Input $(name(input))")
    end
end
"Sets the flag"
Base.setindex!(input::DFInput, dat, key) = setflags!(input, key => dat)

"""
    setflags!(input::DFInput, flags...; print=true)

Sets the specified flags in the input.
"""
function setflags!(input::DFInput{T}, flags...; print=true) where T
    found_keys = Symbol[]
    for (flag, value) in flags
        flag_type = flagtype(input, flag)
        if flag_type !== Nothing
            !(flag in found_keys) && push!(found_keys, flag)
            try
                if isa(value, AbstractVector{<:AbstractVector}) && flag_type <: AbstractVector
                    value = [convert.(eltype(flag_type), v) for v in value]
                else
                    value = convert(flag_type, value)
                end
            catch
                print && ( @warn "Filename '$(name(input))':\n  Could not convert '$value' into '$flag_type'.\n    Flag '$flag' not set.\n")
                continue
            end
            old_data = haskey(input.flags, flag) ? input.flags[flag] : ""
            input.flags[flag] = value
            print && (@info "$(name(input)):\n  -> $flag:\n      $old_data set to: $value\n")
        else
            print && @warn "Flag $flag was ignored since it could not be found in the allowed flags for input $(name(input))."
        end
    end
    return found_keys, input
end

"""
    rmflags!(input::DFInput, flags...)

Remove the specified flags.
"""
function rmflags!(input::DFInput, flags...; print=true)
    for flag in flags
        if haskey(input.flags, flag)
            pop!(input.flags, flag, false)
            print && (@info "Removed flag '$flag' from input '$(name(input))'")
        end
    end
    return input
end

data(input::DFInput, n) = getfirst(x-> name(x) == n, data(input))

"""
    setdata!(input::DFInput, block_name::Symbol, new_block_data; option=nothing, print=true)

sets the data of the specified 'InputData' to the new data. Optionally also sets the 'InputData' option.
"""
function setdata!(input::DFInput, block_name::Symbol, new_block_data; option=nothing, print=true)
    setd = false
    for data_block in input.data
        if data_block.name == block_name
            if typeof(data_block.data) != typeof(new_block_data)
                if print @warn "Overwritten data of type '$(typeof(data_block.data))' with type '$(typeof(new_block_data))'." end
            end
            old_data        = data_block.data
            data_block.data = new_block_data
            data_block.option = option == nothing ? data_block.option : option
            if print
                @info "Block data '$(data_block.name)' in input  '$(name(input))' is now:\n\t$(string(data_block.data)) \n\toption: $(data_block.option)\n"
            end
            setd = true
        end
    end
    if !setd
        setdata!(input, InputData(block_name, option, new_block_data))
        setd = true
    end
    return setd, input
end

"""
    setdataoption!(input::DFInput, name::Symbol, option::Symbol; print=true)

Sets the option of specified data.
"""
function setdataoption!(input::DFInput, name::Symbol, option::Symbol; print=true)
    for data in input.data
        if data.name == name
            old_option  = data.option
            data.option = option
            print && (@info "Option of InputData '$(data.name)' in input '$(name(input))' set from '$old_option' to '$option'")
        end
    end
    return input
end

execs(input::DFInput, exec::String) =
	filter(x -> occursin(exec, x.exec), input.execs)

exec(input::DFInput, exec::String) =
	getfirst(x -> occursin(exec, x.exec), input.execs)

execflags(input::DFInput, exec::String) =
	[x.exec => x.flags for x in execs(input, exec)]

function setexecflags!(input::DFInput, exec::String, flags...)
	for e in execs(input, exec)
		setflags!(e, flags...)
	end
end

setexecdir!(input::DFInput, exec, dir) =
	setexecdir!.(execs(input, exec), dir)

rmexecflags!(input::DFInput, exec::String, flags...) =
	rmflags!.(execs(input, exec), flags...)

runcommand(input::DFInput) = input.execs[1]


"Returns the outputdata for the input."
function outputdata(input::DFInput; print=true, overwrite=true)
    if hasoutput(input)
        if !overwrite && !isempty(outdata(input))
            return outdata(input)
        else
            input.outdata = readoutput(input)
            return input.outdata
        end
    end
    print && (@warn "No output data or output file found for input: $(name(input)).")
    return SymAnyDict()
end
#---------- Extended Interaction ------------#

"""
    setkpoints!(input::DFInput, k_grid)

Sets the kpoints of the input. Will automatically generate the kgrid values if necessary.
"""
function setkpoints!(input::DFInput{Wannier90}, k_grid::NTuple{3, Int}; print=true)
    setflags!(input, :mp_grid => [k_grid...], print=print)
    setdata!(input, :kpoints, kgrid(k_grid..., :wan), print=print)
    return input
end

function setkpoints!(input::DFInput{QE}, k_grid::NTuple{3, Int}; print=true) #nscf

    calc = flag(input, :calculation)
    print && calc != "nscf" && (@warn "Expected calculation to be 'nscf'.\nGot $calc.")
    setdata!(input, :k_points, kgrid(k_grid..., :nscf), option = :crystal, print=print)
    prod(k_grid) > 100 && setflags!(input, :verbosity => "high", print=print)
    return input
end

function setkpoints!(input::DFInput{QE}, k_grid::NTuple{6, Int}; print=true) #scf
    calc = flag(input, :calculation)
    print && calc != "scf" && !occursin("relax", calc) && (@warn "Expected calculation to be scf, vc-relax, relax.\nGot $calc.")
    setdata!(input, :k_points, [k_grid...], option = :automatic, print=print)
    prod(k_grid[1:3]) > 100 && setflags!(input, :verbosity => "high", print=print)
    return input
end

function setkpoints!(input::DFInput{QE}, k_grid::Vector{NTuple{4, T}}; print=true, k_option=:crystal_b) where T<:AbstractFloat
    calc = flag(input, :calculation)
    print && calc != "bands" && (@warn "Expected calculation to be bands, got $calc.")
    @assert in(k_option, [:tpiba_b, :crystal_b, :tpiba_c, :crystal_c]) error("Only $([:tpiba_b, :crystal_b, :tpiba_c, :crystal_c]...) are allowed as a k_option, got $k_option.")
    if k_option in [:tpiba_c, :crystal_c]
        @assert length(k_grid) == 3 error("If $([:tpiba_c, :crystal_c]...) is selected the length of the k_points needs to be 3, got length: $(length(k_grid)).")
    end
    num_k = 0.0
    for k in k_grid
        num_k += k[4]
    end
    if num_k > 100.
        setflags!(input, :verbosity => "high", print=print)
        if print
            @info "Verbosity is set to high because num_kpoints > 100,\n
                       otherwise bands won't get printed."
        end
    end
    setdata!(input, :k_points, k_grid, option=k_option, print=print)
    return input
end

"""
    readbands(input::DFInput)

Tries to read bands from the outputfile of the input.
"""
function readbands(input::DFInput)
    hasoutput_assert(input)
    to = readoutput(input)
    if haskey(to, :bands)
        return to[:bands]
    elseif haskey(to, :bands_up)
        return (up=to[:bands_up], down=to[:bands_down])
    else
        error("No bands found in $(name(input)).")
    end
end

"""
    readfermi(input::DFInput)

Tries to read the fermi level from the outputfile of the input.
"""
function readfermi(input::DFInput)
    hasoutput_assert(input)
    to = readoutput(input)
    if haskey(to, :fermi)
        return to[:fermi]
    else
        error("No fermi found in $(name(input)).")
    end
end

"""
    setwanenergies!(waninput::DFInput{Wannier90}, structure::AbstractStructure, nscf::DFInput , Emin::Real; Epad=5.0)

Automatically calculates and sets the wannier energies. This uses the projections,
`Emin` and the output of the nscf calculation to infer the other limits.
`Epad` allows one to specify the padding around the inner and outer energy windows
"""
function setwanenergies!(input::DFInput{Wannier90}, structure::AbstractStructure, nscf::DFInput, Emin::Real; Epad=5.0)
    hasoutput_assert(nscf)
    iscalc_assert(nscf, "nscf")
    hasprojections_assert(structure)

    bands = readbands(nscf)
    nwann = nprojections(structure)
    @info "num_wann=$nwann (inferred from provided projections)."

    if length(bands) == 2
        num_bands = length(bands[1])
        if hasflag(input, :spin) && input[:spin] == "up"
            winmin, frozmin, frozmax, winmax = wanenergyranges(Emin, nwann, bands.up, Epad)
        elseif hasflag(input, :spin) && input[:spin] == "down"
            winmin, frozmin, frozmax, winmax = wanenergyranges(Emin, nwann, bands.down, Epad)
        end
    else
        num_bands = length(bands)
        winmin, frozmin, frozmax, winmax = wanenergyranges(Emin, nwann, bands, Epad)
    end

    setflags!(input, :dis_win_min => winmin, :dis_froz_min => frozmin, :dis_froz_max => frozmax, :dis_win_max => winmax, :num_wann => nwann, :num_bands=>num_bands;print=false)
    return input
end

#-------- Generating new DFInputs ---------- #

function input_from_kpoints(template::DFInput, newname, kpoints, newflags...)
    newcalc = DFInput(template, newname, newflags...)
    setname!(newcalc, newname)
    setkpoints!(newcalc, kpoints, print=false)
    return newcalc
end

"""
    gencalc_scf(template::DFInput, kpoints::NTuple{6, Int}, newflags...; name="scf")

Uses the information from the template and supplied kpoints to generate an scf input.
Extra flags can be supplied which will be set for the generated input.
"""
gencalc_scf(template::DFInput, kpoints::NTuple{6, Int}, newflags...; name="scf") =
    input_from_kpoints(template, name, kpoints, :calculation => "scf", newflags...)

"""
    gencalc_bands(template::DFInput, kpoints::Vector{NTuple{4}}, newflags...; name="bands")
    gencalc_bands(job::DFJob, kpoints::Vector{NTuple{4}}, newflags...; name="bands", template_name="scf")

Uses the information from the template and supplied kpoints to generate a bands input.
Extra flags can be supplied which will be set for the generated input.
"""
gencalc_bands(template::DFInput, kpoints::Vector{<:NTuple{4}}, newflags...; name="bands") =
    input_from_kpoints(template, name, kpoints, :calculation => "bands", newflags...)


"""
    gencalc_nscf(template::DFInput, kpoints::NTuple{3, Int}, newflags...; name="nscf")
    gencalc_nscf(job::DFJob, kpoints::NTuple{3, Int}, newflags...; name="nscf", template_name="scf")

Uses the information from the template and supplied kpoints to generate an nscf input.
Extra flags can be supplied which will be set for the generated input.
"""
gencalc_nscf(template::DFInput, kpoints::NTuple{3, Int}, newflags...; name="nscf") =
    input_from_kpoints(template, name, kpoints, :calculation => "nscf", newflags...)

"""
    gencalc_projwfc(template::DFInput, Emin, Emax, DeltaE, newflags...; name="projwfc")
    gencalc_projwfc(job::DFJob, Emin, Emax, DeltaE, newflags...; name="projwfc", template_name="nscf")

Uses the information from the template and supplied kpoints to generate a projwfc.x input.
Extra flags can be supplied which will be set for the generated input.
"""
function gencalc_projwfc(template::DFInput, Emin, Emax, DeltaE, extraflags...; name="projwfc")
    occflag = flag(template, :occupations)
    ngauss  = 0
    if occflag == "smearing"
        smearingflag = flag(template, :smearing)
        if smearingflag ∈ ("methfessel-paxton", "m-p", "mp")
            ngauss = 1
        elseif smearingflag ∈ ("marzari-vanderbilt", "cold", "m-v", "mv")
            ngauss = -1
        elseif smearingflag ∈ ("fermi-dirac", "f-d", "fd")
            ngauss = -99
        end
    end
    tdegaussflag = flag(template, :degauss)
    degauss = tdegaussflag != nothing ? tdegaussflag : 0.0
    if length(execs(template)) == 2
        excs = [execs(template)[1], Exec("projwfc.x", execs(template)[end].dir)]
    else
        excs = [Exec("projwfc.x", execs(template)[end].dir)]
    end

    out = DFInput(template, name, excs=excs)
    setname!(out, "projwfc")
    empty!(out.flags)
    setflags!(out, :Emin => Emin, :Emax => Emax, :DeltaE => DeltaE,
              :ngauss => ngauss, :degauss => degauss, print=false)
    setflags!(out, extraflags...)
    return out
end

"""
    gencalc_wan(nscf::DFInput{QE}, structure::AbstractStructure, Emin, wanflags...;
                Epad     = 5.0,
                wanexec  = Exec("wannier90.x", ""))

Generates a Wannier90 input to follow on the supplied `nscf` calculation. It uses the projections defined in the `structure`, and starts counting the required amount of bands from `Emin`.
The `nscf` needs to have a valid output since it will be used in conjunction with `Emin` to find the required amount of bands and energy window for the Wannier90 calculation.
"""
function gencalc_wan(nscf::DFInput{QE}, structure::AbstractStructure, Emin, wanflags...;
                     Epad     = 5.0,
                     wanexec  = Exec("wannier90.x", ""))

    hasoutput_assert(nscf)
    iscalc_assert(nscf, "nscf")
    hasprojections_assert(structure)
    if iscolin(structure)
        wannames = ["wanup", "wandn"]
        @info "Spin polarized calculation found (inferred from nscf input)."
    else
        wannames = ["wan"]
    end

    if flag(nscf, :nosym) != true
        @info "'nosym' flag was not set in the nscf calculation.
                If this was not intended please set it and rerun the nscf calculation.
                This generally gives errors because of omitted kpoints, needed for pw2wannier90.x"
    end
    wanflags = wanflags != nothing ? SymAnyDict(wanflags) : SymAnyDict()

    nwann = nprojections(structure)
    @info "num_wann=$nwann (inferred from provided projections)."
    wanflags[:num_wann]  = nwann
    kpoints = data(nscf, :k_points).data
    wanflags[:mp_grid] = kakbkc(kpoints)
    @info "mp_grid=$(join(wanflags[:mp_grid]," ")) (inferred from nscf input)."
    wanflags[:preprocess] = true
	isnoncolin(structure) && (wanflags[:spinors] = true)
    kdata = InputData(:kpoints, :none, [k[1:3] for k in kpoints])

    waninputs = DFInput{Wannier90}[]
    for  wanfil in wannames
        push!(waninputs, DFInput{Wannier90}(wanfil, dir(nscf), copy(wanflags), [kdata], [Exec(), wanexec], true))
    end

    if length(waninputs) > 1
        setflags!(waninputs[1], :spin => "up")
        setflags!(waninputs[2], :spin => "down")
    end

    map(x -> setwanenergies!(x, structure, nscf, Emin; Epad=Epad), waninputs)
    return waninputs
end

"""
    gencalc_wan(nscf::DFInput{QE}, structure::AbstractStructure, projwfc::DFInput{QE}, threshold::Real, wanflags...; kwargs...)

Generates a wannier calculation, that follows on the `nscf` calculation. Instead of passing Emin manually, the output of a projwfc.x run
can be used together with a `threshold` to determine the minimum energy such that the contribution of the
projections to the DOS is above the `threshold`.
"""
function gencalc_wan(nscf::DFInput{QE}, structure::AbstractStructure, projwfc::DFInput{QE}, threshold::Real, args...; kwargs...)
    hasexec_assert(projwfc, "projwfc.x")
    hasoutput_assert(projwfc)
    Emin = Emin_from_projwfc(structure, projwfc, threshold)
    gencalc_wan(nscf, structure, Emin, args...; kwargs...)
end

"""
    isconverged(input::DFInput{QE})

Returns whether an `scf` calculation was converged.
"""
function isconverged(input::DFInput{QE})
    hasoutput_assert(input)
    iscalc_assert(input, "scf")
    return outputdata(input)[:converged]
end

"""
    setname!(input::DFInput{QE}, name::AbstractString)

Sets the name of `input` to `name` and updates `input.infile` and `input.outfile` to conform
with the new name.
"""
function setname!(input::DFInput{QE}, name::AbstractString; print=true)
    input.name = name
    input.infile = name * splitext(infilename(input))[2]
    input.outfile = name * splitext(outfilename(input))[2]
    print && @info "\ninput.name = $name\ninput.infile = $(infilename(input))\ninput.outfile = $(outfilename(input))"
    name
end
