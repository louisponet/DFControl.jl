#------------- Basic Functionality ---------------#
"""
    save(job::DFJob)

Saves a DFJob, it's job file and all it's input files.
"""
function save(job::DFJob, local_dir=job.local_dir; kwargs...)
    local_dir = local_dir != "" ? local_dir : error("Please specify a valid local_dir!")
    set_localdir!(job, local_dir)
    if isempty(job.name)
        @warn "Job had no name, changed it to: noname"
        job.name = "noname"
    end
    maybe_register_job(job)
    maybe_increment_version(job)
    sanitize_cutoffs!(job)
    sanitize_pseudos!(job)
    sanitize_magnetization!(job)
    sanitize_projections!(job)
    sanitize_flags!(job)
    save_metadata(job)
    return writejobfiles(job; kwargs...)
end

"""
    submit(job::DFJob; server=job.server, server_dir=job.server_dir, rm_prev=true, kwargs...)

Saves the job locally, and then either runs it locally using `qsub` (when `job.server == "localhost"`) or sends it to the specified `job.server` in `job.server_dir`, and submits it using `qsub` on the server. Kwargs get passed through to `save(job; kwargs...)`. If `rm_prev == true` previous `job.tt` output files will be removed.
"""
function submit(job::DFJob; server=job.server, server_dir=job.server_dir, rm_prev=true, kwargs...)
    save(job; kwargs...)
    job.server = server
    try
        job.metadata[:slurmid] = qsub(job; rm_prev=rm_prev)
        save_metadata(job)
    catch
        try
            job.metadata[:slurmid] = sbatch(job; rm_prev=rm_prev)
            save_metadata(job)
        catch
            run(job; rm_prev=rm_prev)
        end
    end
end

"""
    abort(job::DFJob)

Will look for the job id inside it's metadata and try to remove it from the server queue.
If the lastrunning input happened to be a QE input, the correct abort file will be written.
If it's Wannier90 the job will be brutally removed from the slurm queue.
"""
function abort(job::DFJob)
    lastrunning = runninginput(job)
    if lastrunning == nothing
        error("Is this job running?")
    end
    if package(lastrunning) == QE
        length(inputs(job, QE)) > 1 &&
            @warn "It's absolutely impossible to guarantee a graceful abort of a multi job script with QE."

        abortpath = writeabortfile(job, lastrunning)
        while ispath(abortpath)
            continue
        end
        qdel(job)
    else
        qdel(job)
    end
end

"""
    set_flow!(job::DFJob, should_runs...)

Sets whether or not calculations should be run. Calculations are specified using their indices.
"""
function set_flow!(job::DFJob, should_runs...)
    for (name, run) in should_runs
        for input in inputs(job, name)
            input.run = run
        end
    end
    return job
end

"""
    set_headerword!(job::DFJob, old_new::Pair{String, String})


Replaces the specified word in the header with the new word.
"""
function set_headerword!(job::DFJob, old_new::Pair{String, String}; print=true)
    for (i, line) in enumerate(job.header)
        if occursin(first(old_new), line)
            job.header[i] = replace(line, old_new)
            s = """Old line:
            $line
            New line:
            $(job.header[i])
            """
            print && (@info s)
        end
    end
    return job
end

function progressreport(job::DFJob; kwargs...)
    dat = outputdata(job; kwargs...)
    plotdat = SymAnyDict(:fermi=>0.0)
    for (n, d) in dat
        i = input(job, n)
        if isbands(i) && haskey(d, :bands)
            plotdat[:bands] = d[:bands]
        elseif isscf(i) || isnscf(i)
            haskey(d, :fermi) && (plotdat[:fermi] = d[:fermi])
            haskey(d, :accuracy) && (plotdat[:accuracy] = d[:accuracy])
        end
    end
    return plotdat
end

"""
Sets the server dir of the job.
"""
function set_serverdir!(job, dir)
    job.server_dir = dir
    return job
end

"""
    set_localdir!(job::DFJob, dir::String; copy=false)
    
Sets the directory where the job will be saved or ran. If necessary the directory will be created.
If `copy` is set to `true`, all previous inputs and output files of the current job version
(i.e. those in the main job directory) will be copied to the new directory, including the
`outputs` directory with temporary files created during jobs runs.
"""
function set_localdir!(job::DFJob, dir::String; copy=false)
    if !isabspath(dir)
        dir = abspath(dir)
    end
    if !ispath(dir)
        mkpath(dir)
        @info "$dir did not exist, it was created."
    end
    if copy
        cp(job, dir; temp=true)
    end
    job.local_dir = dir
    for i in inputs(job)
        set_dir!(i, dir)
    end
    return job
end


#-------------- Basic Interaction with DFInputs inside the DFJob ---------------#
set_name!(job::DFJob, oldn, newn; kwargs...) = set_name!(input(job, oldn), newn; kwargs...)
Base.insert!(job::DFJob, index::Int, input::DFInput) = insert!(job.inputs, index, input)
Base.push!(job::DFJob, input::DFInput) = push!(job.inputs, input)
Base.pop!(job::DFJob) = pop!(job.inputs)

Base.append!(job::DFJob, args...) = append!(job.inputs, args...)

Base.getindex(job::DFJob, i::Integer) = inputs(job)[i]
Base.length(job::DFJob) = length(inputs(job))
Base.lastindex(job::DFJob) = length(job)

"""Access an input inside the job using it's name. E.g `job["scf"]`"""
function Base.getindex(job::DFJob, id::String)
    tmp = getfirst(x -> name(x)==id, inputs(job))
    if tmp != nothing
        return tmp
    else
        error("No Input with name $id")
    end
end


"""Searches through the inputs for the requested flag.
If a flag was found the input and value of the flag will be added to the returned Dict."""
function Base.getindex(job::DFJob, flg::Symbol)
    outdict = Dict()
    for i in inputs(job)
        tfl = flag(i, flg)
        if tfl != nothing
            outdict[name(i)] = tfl
        end
    end
    return outdict
end

"Set one flag in all the appropriate inputs. E.g `job[:ecutwfc] = 23.0`"
function Base.setindex!(job::DFJob, dat, key::Symbol)
    for input in inputs(job)
        input[key] = dat
    end
end

Base.getindex(job::DFJob, el::Element) = job.structure[el]

"Fuzzily search inputs in the job whose name contain the fuzzy."
searchinputs(job::DFJob, fuzzy::AbstractString) = inputs(job, fuzzy, true)

"Fuzzily search the first input in the job whose name contains the fuzzy."
searchinput(job::DFJob,  fuzzy::AbstractString) = input(job, fuzzy)

"Returns all inputs from a certain package."
searchinputs(job::DFJob, ::Type{P}) where {P <: Package} = filter(x->package(x) == P, inputs(job))


"""
    set_flags!(job::DFJob, inputs::Vector{<:DFInput}, flags...; print=true)

Sets the flags in the names to the flags specified.
This only happens if the specified flags are valid for the names.
If necessary the correct control block will be added to the calculation (e.g. for QEInputs).

The values that are supplied will be checked whether they are valid.
"""
function set_flags!(job::DFJob, inputs::Vector{<:DFInput}, flags...; print=true)
    found_keys = Symbol[]

    for calc in inputs
        t_, = set_flags!(calc, flags..., print=print)
        push!(found_keys, t_...)
    end
    nfound = setdiff([k for (k, v) in flags], found_keys)
    if print && length(nfound) > 0
        f = length(nfound) == 1 ? "flag" : "flags"
        dfprintln("$f '$(join(":" .* String.(setdiff(flagkeys, found_keys)),", "))' were not found in the allowed input variables of the specified inputs!")
    end
    return job
end
set_flags!(job::DFJob, flags...;kwargs...) =
    set_flags!(job, inputs(job), flags...;kwargs...)
set_flags!(job::DFJob, name::String, flags...; fuzzy=true, kwargs...) =
    set_flags!(job, inputs(job, name, fuzzy), flags...; kwargs...)

""" data(job::DFJob, name::String, dataname::Symbol)

Looks through the calculation filenames and returns the data with the specified symbol.
"""
data(job::DFJob, name::String, dataname::Symbol) =
    data(input(job, name), dataname)

"""
    set_data!(job::DFJob, inputs::Vector{<:DFInput}, dataname::Symbol, data; option=nothing)

Looks through the calculation filenames and sets the data of the datablock with `data_block_name` to `new_block_data`.
if option is specified it will set the block option to it.
"""
function set_data!(job::DFJob, inputs::Vector{<:DFInput}, dataname::Symbol, data; kwargs...)
    set_data!.(inputs, dataname, data; kwargs...)
    return job
end
set_data!(job::DFJob, name::String, dataname::Symbol, data; fuzzy=true, kwargs...) =
    set_data!(job, inputs(job, name, fuzzy), dataname, data; kwargs...)

"""
    set_dataoption!(job::DFJob, names::Vector{String}, dataname::Symbol, option::Symbol)

sets the option of specified data in the specified inputs.
"""
function set_dataoption!(job::DFJob, names::Vector{String}, dataname::Symbol, option::Symbol; kwargs...)
    set_dataoption!.(inputs(job, names), dataname, option; kwargs...)
    return job
end
set_dataoption!(job::DFJob, n::String, name::Symbol, option::Symbol; kw...) =
    set_dataoption!(job, [n], name, option; kw...)

"""
    set_dataoption!(job::DFJob, name::Symbol, option::Symbol)

sets the option of specified data block in all calculations that have the block.
"""
set_dataoption!(job::DFJob, n::Symbol, option::Symbol; kw...) =
    set_dataoption!(job, name.(inputs(job)), n, option; kw...)

"""
    rm_flags!(job::DFJob, inputs::Vector{<:DFInput}, flags...)

Looks through the input names and removes the specified flags.
"""
function rm_flags!(job::DFJob, inputs::Vector{<:DFInput}, flags...; kwargs...)
    rm_flags!.(inputs, flags...; kwargs...)
    return job
end
rm_flags!(job::DFJob, name::String, flags...; fuzzy=true, kwargs...) =
    rm_flags!(job, inputs(job, name, fuzzy), flags...; kwargs...)
rm_flags!(job::DFJob, flags...; kwargs...) =
    rm_flags!(job, inputs(job), flags...; kwargs...)

"Returns the executables attached to a given input."
execs(job::DFJob, name) =
    execs(input(job, name))

"""
    set_execflags!(job::DFJob, exec, flags...)

Goes through the calculations of the job and sets the exec flags to the specified ones.
"""
function set_execflags!(job::DFJob, exec, flags...)
	for i in job.inputs
	    set_execflags!(i, exec, flags...)
    end
end

"""
    rmexecflags!(job::DFJob, exec, flags...)

Goes through the calculations of the job and removes the specified `exec` flags.
"""
rmexecflags!(job::DFJob, exec, flags...) =
    rmexecflags!.(job.inputs, (exec, flags)...)

"Sets the directory of the specified executable."
set_execdir!(job::DFJob, exec, dir) =
    set_execdir!.(job.inputs, exec, dir)


"Finds the output files for each of the inputs of a job, and groups all found data into a dictionary."
function outputdata(job::DFJob, inputs::Vector{DFInput}; print=true, onlynew=false)
    datadict = Dict()
    stime = starttime(job)
    for input in inputs
        newout = hasnewout(input, stime)
        if onlynew && !newout
            continue
        end
        tout = outputdata(input; print=print, overwrite=newout)
        if !isempty(tout)
            datadict[name(input)] = tout
        end
    end
    return datadict
end
outputdata(job::DFJob; kwargs...) = outputdata(job, inputs(job); kwargs...)
outputdata(job::DFJob, names::String...; kwargs...) =
    outputdata(job, inputs(job, names); kwargs...)
function outputdata(job::DFJob, n::String; fuzzy=true, kwargs...)
    dat = outputdata(job, inputs(job, n, fuzzy); kwargs...)
    if dat != nothing && haskey(dat, name(input(job, n)))
        return dat[name(input(job, n))]
    end
end

#------------ Specialized Interaction with DFInputs inside DFJob --------------#
"""
    set_kpoints!(job::DFJob, name::String, k_points)

sets the data in the k point `DataBlock` inside the specified inputs.
"""
function set_kpoints!(job::DFJob, n::String, k_points; print=true)
    for calc in inputs(job, n)
        set_kpoints!(calc, k_points, print=print)
    end
    return job
end


"Reads throught the pseudo files and tries to figure out the correct cutoffs"
set_cutoffs!(job::DFJob) = 
    set_cutoffs!.(inputs(job), find_cutoffs(job)...)


"""
    set_wanenergies!(job::DFJob, nscf::DFInput{QE}, Emin::Real; Epad=5.0)

Will set the energy window limits correctly according to the projections specified in the
structure of the job and the specified Emin. The output of `nscf` will be used to determine the
DOS, and what the size of the frozen window needs to be to fit enough bands inside it,
depending on the projections.
"""
function set_wanenergies!(job::DFJob, nscf::DFInput, Emin::Real; Epad=5.0)
    wancalcs = searchinputs(job, Wannier90)
    @assert length(wancalcs) != 0 "Job ($(job.name)) has no Wannier90 calculations, nothing to do."
    map(x->set_wanenergies!(x, structure(job), nscf, Emin; Epad=Epad), wancalcs)
    return job
end

"""
    set_wanenergies!(job::DFJob, nscf::DFInput{QE}, projwfc::DFInput{QE}, threshold::Real; Epad=5.0)

Will set the energy window limits correctly according to the projections specified in the
structure of the job. The output of `projwfc` and the `threshold` will be used to determine
the minimum limit of the frozen energy window such that the interesting DOS of inside it exceeds
the threshold. `nscf` will be used to determine the DOS, and what the upper limit of the frozen window
needs to be to fit enough bands inside it, depending on the projections.
"""
function set_wanenergies!(job::DFJob, nscf::DFInput, projwfc::DFInput, threshold::Real; Epad=5.0)
    hasoutput_assert(projwfc)
    @assert isprojwfc(projwfc) "Please specify a valid projwfc calculation."
    @assert isnscf(nscf) "Please specify a valid nscf calculation."
    Emin = Emin_from_projwfc(job.structure, projwfc, threshold)
    set_wanenergies!(job, nscf, Emin; Epad=Epad)
end

"""
    set_wanenergies!(job::DFJob, min_window_determinator::Real; kwargs...)

Sets the energy windows of wannier calculations based on the `job`.
When a projwfc calculation is present in the `job`, `min_window_determinator` will be used to
determine the threshold value for including a band in the window based on the projections, otherwise
it will be used as the `Emin` value from which to start counting the number of bands needed for all
projections.
"""
function set_wanenergies!(job::DFJob, min_window_determinator::Real; kwargs...)
    nscf_input = getfirst(isnscf, inputs(job))
    projwfc_input = getfirst(isprojwfc, inputs(job))
    if projwfc_input === nothing || !hasoutput(projwfc_input)
        @info "No projwfc input found with valid output, using $min_window_determinator as Emin"
        return set_wanenergies!(job, nscf_input, min_window_determinator; kwargs...)
    else
        @info "Valid projwfc output found, using $min_window_determinator as the dos threshold."
        return set_wanenergies!(job, nscf_input, projwfc_input, min_window_determinator; kwargs...)
    end
end

#--------------- Interacting with the Structure inside the DFJob ---------------#
"Returns the ith atom with id `atsym`."
atom(job::DFJob, atsym::Symbol, i=1) = filter(x -> x.name == atsym, atoms(job))[i]

"""
    atoms(job::DFJob)

Returns the atoms inside the structure of the job.
"""
atoms(job::DFJob) = atoms(job.structure)
atoms(job::DFJob, args...) = atoms(job.structure, args...)
atoms(f::Function, job::DFJob) = atoms(f, job.structure)

"""job.structure.atoms = atoms"""
set_atoms!(job::DFJob, atoms::Vector{<:AbstractAtom}) =
    job.structure.atoms = atoms

#automatically sets the cell parameters for the entire job, implement others
"""
    set_cell!(job::DFJob, cell_::Mat3)

sets the cell parameters of the structure in the job.
"""
function set_cell!(job::DFJob, cell_::Mat3)
    job.structure.cell = cell_
    return job
end

"sets the pseudopotentials to the specified one in the default pseudoset."
set_pseudos!(job::DFJob, set::Symbol, specifier::String=""; kwargs...) = 
    set_pseudos!(job.structure, set, specifier; kwargs...)

"sets the pseudopotentials for the atom with name `atsym` to the specified one in the default pseudoset."
set_pseudos!(job::DFJob, atsym::Symbol, set::Symbol, specifier::String=""; kwargs...) =
	set_pseudos!(job.structure, atsym, set, specifier; kwargs...)

"sets the pseudopotentials to the specified one in the default pseudoset."
set_pseudos!(job::DFJob, at_pseudos::Pair{Symbol, Pseudo}...; kwargs...) = 
    set_pseudos!(job.structure, at_pseudos...; kwargs...)

"Returns the projections inside the job for the specified `i`th atom in the job with id `atsym`."
projections(job::DFJob, atsym::Symbol, i=1) = projections(atom(job, atsym, i))

"Returns all the projections inside the job."
projections(job::DFJob) = projections(structure(job))

"""
sets the projections of the specified atoms inside the job structure.
"""
function set_projections!(job::DFJob, projections...; kwargs...)
    socid = findfirst(issoc, inputs(job))
    set_projections!(job.structure, projections...; soc=socid !== nothing, kwargs...)
end

for hub_param in (:U, :J0, :α, :β)
	f = Symbol("set_Hubbard_$(hub_param)!")
	str = "$hub_param"
	@eval begin
		"""
			$($(f))(job::DFJob, ats_$($(str))s::Pair{Symbol, <:AbstractFloat}...; print=true)

		Set the Hubbard $($(str)) parameter for the specified atoms.

		Example:
			`$($(f))(job, :Ir => 2.1, :Ni => 1.0, :O => 0.0)`
		"""
		function $f(job::DFJob, $(hub_param)::Pair{Symbol, <:AbstractFloat}...; print=true)
			for (atsym, val) in $(hub_param)
				$f.(atoms(job, atsym), val; print=print)
			end
		end
		export $f
	end
end

"""
	set_Hubbard_J!(job::DFJob, ats_Js::Pair{Symbol, Vector{<:AbstractFloat}}...; print=true)

Set the Hubbard J parameter for the specified atom.

Example:
	`set_Hubbard_J(job, :Ir => [2.1], :Ni => [1.0])'
"""
function set_Hubbard_J!(job::DFJob, ats_Js::Pair{Symbol, <:Vector{<:AbstractFloat}}...; print=true)
	for (atsym, val) in ats_Js
		set_Hubbard_J!.(atoms(job, atsym), (val,); print=print)
	end
end

export set_Hubbard_J!


"Rescales the unit cell."
scale_cell!(job::DFJob, s) =
	scale_cell!(job.structure, s)

set_magnetization!(job::DFJob, args...) =
	set_magnetization!(job.structure, args...)

Base.joinpath(job::DFJob, p) =
	joinpath(job.local_dir, p)

create_supercell(job::DFJob, args...) =
	create_supercell(structure(job), args...)

volume(job::DFJob) = volume(structure(job))

"""
	bandgap(job::DFJob, fermi=nothing)

Calculates the bandgap (possibly indirect) around the fermi level.
Uses the first found bands calculation, if there is none it uses the first found nscf calculation.
"""
function bandgap(job::DFJob, fermi=nothing)
	band_calcs = getfirst.([isbands, isnscf, isscf], (inputs(job),))
	if all(x -> x === nothing, band_calcs)
		error("No valid calculation found to calculate the bandgap.\nMake sure the job has either a valid bands or nscf calculation.")
	end
	if fermi === nothing
    	fermi_calcs = getfirst.([isnscf, isscf], (inputs(job),))
    	if all(x -> x === nothing, band_calcs)
    		error("No valid calculation found to extract the fermi level.\nPlease supply the fermi level manually.")
		end
		fermi = maximum(readfermi.(filter(x -> x!==nothing, fermi_calcs)))
	end

	bands = readbands.(filter(x -> x!==nothing, band_calcs))
	return minimum(bandgap.(bands, fermi))
end

"Searches for the first scf or nscf calculation with output, and reads the fermi level from it."
function readfermi(job::DFJob)
    ins = filter(x -> (isscf(x) || isnscf(x)) && hasoutfile(x), inputs(job))
    @assert isempty(ins) !== nothing "Job does not have a valid scf or nscf output."
    for i in ins
        o = outputdata(i)
        if haskey(o, :fermi)
            return o[:fermi]
        end
    end
    @warn "No output files with fermi level found." 
    return 0.0
end

"Searches for the first bands calculation with output, and reads the fermi level from it."
function readbands(job::DFJob)
    input = getfirst(x -> isbands(x) && hasoutfile(x), inputs(job))
    if input === nothing
        input = getfirst(x -> isnscf(x) && hasoutfile(x), inputs(job))
        if input === nothing
            @warn "Job does not have a valid bands output."
            return nothing
        end
        @warn "No bands calculation found, return bands from nscf calculation." 
        return readbands(input)
    end
    return readbands(input)
end

"""
    symmetry_operators(j::DFJob; maxsize=52, tolerance=$DEFAULT_TOLERANCE)
    symmetry_operators(s::Structure; maxsize=52, tolerance=$DEFAULT_TOLERANCE)

Finds and returns all the rotations and translations that are symmetry operators of the structure.
"""
symmetry_operators(j::DFJob; kwargs...) = symmetry_operators(j.structure; kwargs...)

"""
    international(j::DFJob; tolerance=$DEFAULT_TOLERANCE)
    international(s::Structure; tolerance=$DEFAULT_TOLERANCE)

Returns the international symbol of the space group of the structure.
"""
international(j::DFJob; kwargs...) = international(j.structure; kwargs...)

"""
    niggli_reduce(j::DFJob; tolerance=$DEFAULT_TOLERANCE)
    niggli_reduce(s::Structure; tolerance=$DEFAULT_TOLERANCE)

Returns the niggli reduced lattice cell.
"""
niggli_reduce(j::DFJob; kwargs...) = niggli_reduce(j.structure; kwargs...)

for (calc, tn) in zip((:gencalc_bands, :gencalc_nscf, :gencalc_projwfc), ("scf", "scf", "nscf"))
    @eval function $calc(job::DFJob, args...;template_name::String=$tn, kwargs...)
            template = input(job, template_name)
            @assert template !== nothing "No valid input with template_name $template_name found in job."
            return $calc(template, args...; kwargs...)
        end
end

"""
    gencalc_wan(job::DFJob, min_window_determinator::Real, extra_wan_flags...; kwargs...)

Automates the generation of wannier calculations based on the `job`.
When a projwfc calculation is present in the `job`, `min_window_determinator` will be used to
determine the threshold value for including a band in the window based on the projections, otherwise
it will be used as the `Emin` value from which to start counting the number of bands needed for all
projections.
`extra_wan_flags` can be any extra flags for the Wannier90 input such as `write_hr` etc.
"""
function gencalc_wan(job::DFJob, min_window_determinator::Real, extra_wan_flags...; kwargs...)
    nscf_input = getfirst(x -> isnscf(x), inputs(job))
    projwfc_input = getfirst(x -> isprojwfc(x), inputs(job))
    if projwfc_input === nothing || !hasoutput(projwfc_input)
        @info "No projwfc input found with valid output, using $min_window_determinator as Emin"
        return gencalc_wan(nscf_input, job.structure, min_window_determinator, extra_wan_flags...; kwargs...)
    else
        @info "Valid projwfc output found, using $min_window_determinator as the dos threshold."
        return gencalc_wan(nscf_input, job.structure, projwfc_input, min_window_determinator, extra_wan_flags...; kwargs...)
    end
end

#TODO: only for QE 
"Reads the pdos for a particular atom. Only works for QE."  
function pdos(job::DFJob, atsym::Symbol, filter_word="") 
    projwfc = getfirst(isprojwfc, inputs(job)) 
    ats = atoms(job, atsym)
    @assert length(ats) > 0 "No atoms found with name $atsym."
    scf = getfirst(isscf, inputs(job))
    magnetic = any(ismagnetic, atoms(job)) || ismagnetic(scf) 
    soc = issoc(scf)
    return pdos(projwfc, atsym, magnetic, soc, filter_word)
end

pdos(job::DFJob, atom::AbstractAtom, args...) =
    pdos(job, name(atom), args...)

function pdos(job::DFJob, atoms::Vector{AbstractAtom} = atoms(job), args...)
    t_energies, t_pdos = pdos(job, atoms[1], args...)
    for i in 2:length(atoms)
        t1, t2 = pdos(job, atoms[i], args...)
        t_pdos .+= t2
    end
    return (energies=t_energies, pdos=t_pdos)
end

update_geometry!(job::DFJob, new_str::AbstractStructure) =
    update_geometry!(job.structure, new_str)

high_symmetry_kpath(j::DFJob, args...;kwargs...) =
    high_symmetry_kpath(j.structure, args...; kwargs...)

