#------------- Basic Functionality ---------------#
"""
    save(job::DFJob)

Saves a DFJob, it's job file and all it's input files.
"""
function save(job::DFJob, local_dir=job.local_dir; kwargs...)
    local_dir = local_dir != "" ? local_dir : error("Please specify a valid local_dir!")
    if !ispath(local_dir)
        mkpath(local_dir)
        @info "$local_dir did not exist, it was created."
    end
    sanitizeflags!(job)
    job.local_dir = local_dir
    return writejobfiles(job; kwargs...)
end

"""
    submit(job::DFJob; server=job.server, server_dir=job.server_dir, rm_prev=true, kwargs...)

Saves the job locally, and then either runs it locally using `qsub` (when `job.server == "localhost"`) or sends it to the specified `job.server` in `job.server_dir`, and submits it using `qsub` on the server. Kwargs get passed through to `save(job; kwargs...)`. If `rm_prev == true` previous `job.tt` output files will be removed.
"""
function submit(job::DFJob; server=job.server, server_dir=job.server_dir, rm_prev=true, kwargs...)
    save(job; kwargs...)
    job.server = server
    job.metadata[:slurmid] = qsub(job; rm_prev=rm_prev)
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
    setflow!(job::DFJob, should_runs...)

Sets whether or not calculations should be run. Calculations are specified using their indices.
"""
function setflow!(job::DFJob, should_runs...)
    for (name, run) in should_runs
        for input in inputs(job, name)
            input.run = run
        end
    end
    return job
end

"""
    setheaderword!(job::DFJob, old_new::Pair{String, String})


Replaces the specified word in the header with the new word.
"""
function setheaderword!(job::DFJob, old_new::Pair{String, String}; print=true)
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
        if isbandscalc(i) && haskey(d, :bands)
            plotdat[:bands] = d[:bands]
        elseif isscfcalc(i) || isnscfcalc(i)
            haskey(d, :fermi) && (plotdat[:fermi] = d[:fermi])
            haskey(d, :accuracy) && (plotdat[:accuracy] = d[:accuracy])
        end
    end
    return plotdat
end

"""
Sets the server dir of the job.
"""
function setserverdir!(job, dir)
    job.server_dir = dir
    return job
end

"""
Sets the local dir of the job.
"""
function setlocaldir!(job, dir)
    if !isabspath(dir)
        dir = abspath(dir)
    end
    job.local_dir = dir
    for i in inputs(job)
        setdir!(i, dir)
    end
    return job
end


#-------------- Basic Interaction with DFInputs inside the DFJob ---------------#
setname!(job::DFJob, oldn, newn) = (input(job, oldn).name = newn)
Base.insert!(job::DFJob, index::Int, input::DFInput) = insert!(job.inputs, index, input)
Base.push!(job::DFJob, input::DFInput) = push!(job.inputs, input)
Base.pop!(job::DFJob) = pop!(job.inputs)

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

"Fuzzily search inputs in the job whose name contain the fuzzy."
searchinputs(job::DFJob, fuzzy::AbstractString) = inputs(job, fuzzy, true)

"Fuzzily search the first input in the job whose name contains the fuzzy."
searchinput(job::DFJob,  fuzzy::AbstractString) = input(job, fuzzy)

"Returns all inputs from a certain package."
searchinputs(job::DFJob, ::Type{P}) where {P <: Package} = filter(x->package(x) == P, inputs(job))


"""
    setflags!(job::DFJob, inputs::Vector{<:DFInput}, flags...; print=true)

Sets the flags in the names to the flags specified.
This only happens if the specified flags are valid for the names.
If necessary the correct control block will be added to the calculation (e.g. for QEInputs).

The values that are supplied will be checked whether they are valid.
"""
function setflags!(job::DFJob, inputs::Vector{<:DFInput}, flags...; print=true)
    found_keys = Symbol[]

    for calc in inputs
        t_, = setflags!(calc, flags..., print=print)
        push!(found_keys, t_...)
    end
    nfound = setdiff([k for (k, v) in flags], found_keys)
    if print && length(nfound) > 0
        f = length(nfound) == 1 ? "flag" : "flags"
        dfprintln("$f '$(join(":" .* String.(setdiff(flagkeys, found_keys)),", "))' were not found in the allowed input variables of the specified inputs!")
    end
    return job
end
setflags!(job::DFJob, flags...;kwargs...) =
    setflags!(job, inputs(job), flags...;kwargs...)
setflags!(job::DFJob, name::String, flags...; fuzzy=true, kwargs...) =
    setflags!(job, inputs(job, name, fuzzy), flags...; kwargs...)

""" data(job::DFJob, name::String, dataname::Symbol)

Looks through the calculation filenames and returns the data with the specified symbol.
"""
data(job::DFJob, name::String, dataname::Symbol) =
    data(input(job, name), dataname)

"""
    setdata!(job::DFJob, inputs::Vector{<:DFInput}, dataname::Symbol, data; option=nothing)

Looks through the calculation filenames and sets the data of the datablock with `data_block_name` to `new_block_data`.
if option is specified it will set the block option to it.
"""
function setdata!(job::DFJob, inputs::Vector{<:DFInput}, dataname::Symbol, data; kwargs...)
    setdata!.(inputs, dataname, data; kwargs...)
    return job
end
setdata!(job::DFJob, name::String, dataname::Symbol, data; fuzzy=true, kwargs...) =
    setdata!(job, inputs(job, name, fuzzy), dataname, data; kwargs...)

"""
    setdataoption!(job::DFJob, names::Vector{String}, dataname::Symbol, option::Symbol)

sets the option of specified data in the specified inputs.
"""
function setdataoption!(job::DFJob, names::Vector{String}, dataname::Symbol, option::Symbol; kwargs...)
    setdataoption!.(inputs(job, names), dataname, option; kwargs...)
    return job
end
setdataoption!(job::DFJob, n::String, name::Symbol, option::Symbol; kw...) =
    setdataoption!(job, [n], name, option; kw...)

"""
    setdataoption!(job::DFJob, name::Symbol, option::Symbol)

sets the option of specified data block in all calculations that have the block.
"""
setdataoption!(job::DFJob, n::Symbol, option::Symbol; kw...) =
    setdataoption!(job, name.(inputs(job)), n, option; kw...)

"""
    rmflags!(job::DFJob, inputs::Vector{<:DFInput}, flags...)

Looks through the input names and removes the specified flags.
"""
function rmflags!(job::DFJob, inputs::Vector{<:DFInput}, flags...; kwargs...)
    rmflags!.(inputs, flags...; kwargs...)
    return job
end
rmflags!(job::DFJob, name::String, flags...; fuzzy=true, kwargs...) =
    rmflags!(job, inputs(job, name, fuzzy), flags...; kwargs...)
rmflags!(job::DFJob, flags...; kwargs...) =
    rmflags!(job, inputs(job), flags...; kwargs...)

"Returns the executables attached to a given input."
execs(job::DFJob, name) =
    execs(input(job, name))

"""
    setexecflags!(job::DFJob, exec, flags...)

Goes through the calculations of the job and sets the exec flags to the specified ones.
"""
setexecflags!(job::DFJob, exec, flags...) =
    setexecflags!.(job.inputs, (exec, flags)...)
"""
    rmexecflags!(job::DFJob, exec, flags...)

Goes through the calculations of the job and removes the specified `exec` flags.
"""
rmexecflags!(job::DFJob, exec, flags...) =
    rmexecflags!.(job.inputs, (exec, flags)...)

"Sets the directory of the specified executable."
setexecdir!(job::DFJob, exec, dir) =
    setexecdir!.(job.inputs, exec, dir)


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
    datadict
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
    setkpoints!(job::DFJob, name::String, k_points)

sets the data in the k point `DataBlock` inside the specified inputs.
"""
function setkpoints!(job::DFJob, name::String, k_points; print=true)
    for calc in inputs(job, n)
        setkpoints!(calc, k_points, print=print)
    end
    return job
end


"Reads throught the pseudo files and tries to figure out the correct cutoffs"
function setcutoffs!(job::DFJob)
    @assert job.server == "localhost" "Cutoffs can only be automatically set if the pseudo files live on the local machine."
    pseudofiles = map(x->x.name, filter(!isempty, [pseudo(at) for at in atoms(job)]))
    pseudodirs  = map(x->x.dir, filter(!isempty, [pseudo(at) for at in atoms(job)]))
    @assert !isempty(pseudofiles) "No atoms with pseudo files found."
    @assert !isempty(pseudodirs) "No valid pseudo directories found in the inputs."
    maxecutwfc = 0.0
    maxecutrho = 0.0
    for d in pseudodirs
        for f in pseudofiles
            pth = joinpath(d, f)
            if ispath(pth)
                ecutwfc, ecutrho = read_cutoffs_from_pseudofile(pth)
                if ecutwfc != nothing && ecutrho != nothing
                    maxecutwfc = ecutwfc > maxecutwfc ? ecutwfc : maxecutwfc
                    maxecutrho = ecutrho > maxecutrho ? ecutrho : maxecutrho
                end
            end
        end
    end
    setcutoffs!.(inputs(job), maxecutwfc, maxecutrho)
end

"""
    setwanenergies!(job::DFJob, nscf::DFInput{QE}, Emin::Real; Epad=5.0)

Will set the energy window limits correctly according to the projections specified in the
structure of the job and the specified Emin. The output of `nscf` will be used to determine the
DOS, and what the size of the frozen window needs to be to fit enough bands inside it,
depending on the projections.
"""
function setwanenergies!(job::DFJob, nscf::DFInput{QE}, Emin::Real; Epad=5.0)
    hasoutput_assert(nscf)
    iscalc_assert(nscf, "nscf")
    hasprojections_assert(structure(job))

    bands = readbands(nscf)
    wancalcs = searchinputs(job, Wannier90)
    @assert length(wancalcs) != 0 "Job ($(job.name)) has no Wannier90 calculations, nothing to do."
    nbnd = isnoncolincalc(nscf) ? 2 * nprojections(structure(job)) : nprojections(structure(job))
    @info "num_bands=$nbnd (inferred from provided projections)."
    winmin, frozmin, frozmax, winmax = wanenergyranges(Emin, nbnd, bands, Epad)
    map(x->setflags!(x, :dis_win_min => winmin, :dis_froz_min => frozmin, :dis_froz_max => frozmax, :dis_win_max => winmax, :num_wann => nbnd, :num_bands=>length(bands)), wancalcs)
    return job
end

"""
    setwanenergies!(job::DFJob, nscf::DFInput{QE}, projwfc::DFInput{QE}, threshold::Real; Epad=5.0)

Will set the energy window limits correctly according to the projections specified in the
structure of the job. The output of `projwfc` and the `threshold` will be used to determine
the minimum limit of the frozen energy window such that the interesting DOS of inside it exceeds
the threshold. `nscf` will be used to determine the DOS, and what the upper limit of the frozen window
needs to be to fit enough bands inside it, depending on the projections.
"""
function setwanenergies!(job::DFJob, nscf::DFInput{QE}, projwfc::DFInput{QE}, threshold::Real; Epad=5.0)
    hasoutput_assert(projwfc)
    hasexec_assert(projwfc, "projwfc.x")
    Emin = Emin_from_projwfc(job.structure, projwfc, threshold)
    setwanenergies!(job, nscf, Emin; Epad=Epad)
end

#--------------- Interacting with the Structure inside the DFJob ---------------#
"Returns the ith atom with id `atsym`."
atom(job::DFJob, atsym::Symbol, i=1) = filter(x -> x.name == atsym, atoms(job))[i]

"""
    atoms(job::DFJob)

Returns a list the atoms in the structure.
"""
atoms(job::DFJob) = atoms(job.structure)

"""
    setatoms!(job::DFJob, atoms::Dict{Symbol,<:Array{<:Point3,1}}, pseudo_setname=nothing, pseudospecifier=nothing, option=:angstrom)

Sets the data data with atomic positions to the new one. This is done for all calculations in the job that have that data.
If default pseudopotentials are defined, a set can be specified, together with a fuzzy that distinguishes between the possible multiple pseudo strings in the pseudo set.
These pseudospotentials are then set in all the calculations that need it.
All flags which specify the number of atoms inside the calculation also gets set to the correct value.
"""
function setatoms!(job::DFJob, atoms::Vector{<:AbstractAtom}; pseudoset=nothing, pseudospecifier="")
    job.structure.atoms = atoms
    pseudoset!=nothing && setpseudos!(job, pseudoset, pseudospecifier)
    return job
end

#automatically sets the cell parameters for the entire job, implement others
"""
    setcell!(job::DFJob, cell_::Mat3)

sets the cell parameters of the structure in the job.
"""
function setcell!(job::DFJob, cell_::Mat3)
    job.structure.cell = cell_
    return job
end

"sets the pseudopotentials to the specified one in the default pseudoset."
setpseudos!(job::DFJob, set::Symbol, specifier::String=""; kwargs...) = 
    setpseudos!(job.structure, set, specifier; kwargs...)

"sets the pseudopotentials for the atom with name `atsym` to the specified one in the default pseudoset."
setpseudos!(job::DFJob, atsym::Symbol, set::Symbol, specifier::String=""; kwargs...) =
	setpseudos!(job.structure, atsym, set, specifier; kwargs...)

"sets the pseudopotentials to the specified one in the default pseudoset."
setpseudos!(job::DFJob, at_pseudos::Pair{Symbol, Pseudo}...; kwargs...) = 
    setpseudos!(job.structure, at_pseudos...; kwargs...)

"Returns the projections inside the job for the specified `i`th atom in the job with id `atsym`."
projections(job::DFJob, atsym::Symbol, i=1) = projections(atom(job, atsym, i))

"Returns all the projections inside the job."
projections(job::DFJob) = projections(structure(job))

"""
sets the projections of the specified atoms inside the job structure.
"""
setprojections!(job::DFJob, projections...) =
    setprojections!(job.structure, projections...)

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

atoms(job::DFJob, atsym::Symbol) = atoms(job.structure, atsym)

"Rescales the unit cell."
scale_cell!(job::DFJob, s) =
	scale_cell!(job.structure, s)

set_magnetization!(job::DFJob, args...) =
	set_magnetization!(job.structure, args...)

Base.joinpath(job::DFJob, p) =
	joinpath(job.local_dir, p)

create_supercell(job::DFJob, args...) =
	create_supercell(structure(job), args...)
