#------------- Basic Functionality ---------------#
"""
    save(job::DFJob)

Saves the job's calculations and `job.tt` submission script in `job.local_dir`.
Some sanity checks will be performed on the validity of flags, execs, pseudopotentials, etc.
The job will also be registered for easy retrieval at a later stage.

If a previous job is present in the job directory (indicated by a valid job script),
it will be copied to the `.versions` sub directory as the previous version of `job`,
and the version of `job` will be incremented. 
"""
function save(job::DFJob; kwargs...)
    # First we check whether the job is trying to be saved in a archived directory, absolutely not allowed
    @assert !isarchived(job)
        "Not allowed to save a job in a archived directory, please specify a different directory with `set_localdir!"

    # Here we find the main directory, needed for if a job's local dir is a .versions one
    local_dir = main_job_dir(job)
    if ispath(joinpath(local_dir, "job.tt"))
        tj = DFJob(local_dir)
        @assert !isrunning(tj) "Can't save a job in a directory where another is running."
        # We want to first make sure that the previous job in the main directory is safely stored
        cp(tj, joinpath(tj, VERSION_DIR_NAME, "$(tj.version)"); force = true)
    end
    if local_dir != job.local_dir
        # We know for sure it was a previously saved job
        # Now that we have safely stored it we can clean out the directory to then fill
        # it with the files from the job.version
        clean_local_dir!(local_dir)
        cp(job, local_dir; force = true)
    end
    if isempty(job.name)
        @warn "Job had no name, changed it to: noname"
        job.name = "noname"
    end
    
    set_localdir!(job, local_dir) # Needs to be done so the inputs `dir` also changes.
    verify_or_create(local_dir)
    maybe_register_job(job)

    curver = job.version
    job.version = last_version(job) + 1
    @info "Job version: $(curver) => $(job.version)."
    sanitize_cutoffs!(job)
    sanitize_pseudos!(job)
    sanitize_magnetization!(job)
    sanitize_projections!(job)
    sanitize_flags!(job)
    timestamp!(job, now())
    save_metadata(job)
    return writejobfiles(job; kwargs...)
end

"""
    submit(job::DFJob; kwargs...)

First saves the job, then tries to submit the job script through `sbatch job.tt` or `bash job.tt` if the former command fails.
`kwargs...` get passed to `save(job; kwargs...)`.
"""
function submit(job::DFJob; server = job.server, server_dir = job.server_dir, kwargs...)
    save(job; kwargs...)
    job.server = server
    try
        job.metadata[:slurmid] = qsub(job)
        save_metadata(job)
    catch
        try
            job.metadata[:slurmid] = sbatch(job)
            save_metadata(job)
        catch
            pop!(job.metadata, :slurmid, nothing)
            run(job)
        end
    end
end

"""
    abort(job::DFJob)

Will try to remove the job from the scheduler's queue.
If the last running calculation happened to be a `DFCalculation{QE}`, the correct abort file will be written.
For other codes the process is not smooth, and restarting is not guaranteed.
"""
function abort(job::DFJob)
    lastrunning = last_running_calculation(job)
    if lastrunning == nothing
        error("Is this job running?")
    end
    if package(lastrunning) == QE
        length(calculations(job, QE)) > 1 &&
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
    set_flow!(job::DFJob, should_runs::Pair{String, Bool}...)

Sets whether or not calculations should be scheduled to run.
The `name` of each calculation in the job will be checked against the string in each pair of `should_runs`, and the
`calculation.run` will be set accordingly.

Example:
```julia
set_flow!(job, "" => false, "scf" => true)
```
would un-schedule all calculations in the job, and schedule the "scf" and "nscf" calculations to run.
"""
function set_flow!(job::DFJob, should_runs...)
    for (name, run) in should_runs
        for calculation in calculations(job, name)
            calculation.run = run
        end
    end
    return job
end

"""
    set_headerword!(job::DFJob, old_new::Pair{String, String})

Replaces the specified word in the header with the new word.
"""
function set_headerword!(job::DFJob, old_new::Pair{String,String}; print = true)
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
    plotdat = SymAnyDict(:fermi => 0.0)
    for (n, d) in dat
        i = calculation(job, n)
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
    set_localdir!(job::DFJob, dir::AbstractString; copy=false)
    
Sets `job.local_dir` to `dir`. If necessary the directory will be created.
If `copy` is set to `true`, all previous calculations and output files of the current job version
(i.e. those in the main job directory) will be copied to the new directory, including the
`outputs` directory with temporary files created during jobs runs.
"""
function set_localdir!(job::DFJob, dir::AbstractString; copy = false)
    if !isabspath(dir)
        dir = abspath(dir)
    end
    if dir[end] == '/'
        dir = dir[1:end-1]
    end
    if copy
        verify_or_create(dir)
        cp(job, dir; temp = true)
    end
    job.local_dir = dir
    for i in calculations(job)
        set_dir!(i, dir)
    end
    return job
end

#-------------- Basic Interaction with DFCalculations inside the DFJob ---------------#
"""
    insert!(job::DFJob, i::Int, calculation::DFCalculation) = insert!(job.calculations, i, calculation)
"""
function Base.insert!(job::DFJob, index::Int, calculation::DFCalculation)
    return insert!(job.calculations, index, calculation)
end

"""
    push!(job::DFJob, calculation::DFCalculation) = push!(job.calculations, calculation)
"""
Base.push!(job::DFJob, calculation::DFCalculation) = push!(job.calculations, calculation)

"""
    pop!(job::DFJob) = pop!(job.calculations)
"""
Base.pop!(job::DFJob) = pop!(job.calculations)

"""
    append!(job::DFJob, args...) = append!(job.calculations, args...)
"""
Base.append!(job::DFJob, args...) = append!(job.calculations, args...)

Base.getindex(job::DFJob, i::Integer) = calculations(job)[i]
Base.length(job::DFJob) = length(calculations(job))
Base.lastindex(job::DFJob) = length(job)

"""
    getindex(job::DFJob, name::String)
    
Returns the `DFCalculation` with the specified `name`.

    getindex(job::DFJob, i::Integer)
    
Returns the i'th `DFCalculation` in the job.
"""
function Base.getindex(job::DFJob, id::String)
    tmp = getfirst(x -> name(x) == id, calculations(job))
    if tmp != nothing
        return tmp
    else
        error("No Input with name $id")
    end
end

Base.getindex(job::DFJob, el::Element) = job.structure[el]

"""
    set_data!(job::DFJob, calculations::Vector{<:DFCalculation}, dataname::Symbol, data; option=nothing)

Looks through the calculation filenames and sets the data of the datablock with `data_block_name` to `new_block_data`.
if option is specified it will set the block option to it.
"""
function set_data!(job::DFJob, calculations::Vector{<:DFCalculation}, dataname::Symbol,
                   data; kwargs...)
    set_data!.(calculations, dataname, data; kwargs...)
    return job
end
function set_data!(job::DFJob, name::String, dataname::Symbol, data; fuzzy = true,
                   kwargs...)
    return set_data!(job, calculations(job, name, fuzzy), dataname, data; kwargs...)
end

"""
    set_data_option!(job::DFJob, names::Vector{String}, dataname::Symbol, option::Symbol)

sets the option of specified data in the specified calculations.
"""
function set_data_option!(job::DFJob, names::Vector{String}, dataname::Symbol,
                          option::Symbol; kwargs...)
    set_data_option!.(calculations(job, names), dataname, option; kwargs...)
    return job
end
function set_data_option!(job::DFJob, n::String, name::Symbol, option::Symbol; kw...)
    return set_data_option!(job, [n], name, option; kw...)
end

"""
    set_data_option!(job::DFJob, name::Symbol, option::Symbol)

sets the option of specified data block in all calculations that have the block.
"""
function set_data_option!(job::DFJob, n::Symbol, option::Symbol; kw...)
    return set_data_option!(job, name.(calculations(job)), n, option; kw...)
end

"Finds the output files for each of the calculations of a job, and groups all found data into a dictionary."
function outputdata(job::DFJob, calculations::Vector{DFCalculation}; print = true,
                    onlynew = false)
    datadict = Dict()
    stime = starttime(job)
    for calculation in calculations
        newout = hasnewout(calculation, stime)
        if onlynew && !newout
            continue
        end
        tout = outputdata(calculation; print = print, overwrite = newout)
        if !isempty(tout)
            datadict[name(calculation)] = tout
        end
    end
    return datadict
end
outputdata(job::DFJob; kwargs...) = outputdata(job, calculations(job); kwargs...)
function outputdata(job::DFJob, names::String...; kwargs...)
    return outputdata(job, calculations(job, names); kwargs...)
end
function outputdata(job::DFJob, n::String; fuzzy = true, kwargs...)
    dat = outputdata(job, calculations(job, n, fuzzy); kwargs...)
    if dat != nothing && haskey(dat, name(calculation(job, n)))
        return dat[name(calculation(job, n))]
    end
end

#------------ Specialized Interaction with DFCalculations inside DFJob --------------#
"Reads throught the pseudo files and tries to figure out the correct cutoffs"
set_cutoffs!(job::DFJob) = set_cutoffs!.(calculations(job), find_cutoffs(job)...)

"""
    set_wanenergies!(job::DFJob, nscf::DFCalculation{QE}, Emin::Real; Epad=5.0)

Will set the energy window limits correctly according to the projections specified in the
structure of the job and the specified Emin. The output of `nscf` will be used to determine the
DOS, and what the size of the frozen window needs to be to fit enough bands inside it,
depending on the projections.
"""
function set_wanenergies!(job::DFJob, nscf::DFCalculation, Emin::Real; Epad = 5.0)
    wancalcs = calculations(job, Wannier90)
    @assert length(wancalcs) != 0 "Job ($(job.name)) has no Wannier90 calculations, nothing to do."
    map(x -> set_wanenergies!(x, structure(job), nscf, Emin; Epad = Epad), wancalcs)
    return job
end

"""
    set_wanenergies!(job::DFJob, nscf::DFCalculation{QE}, projwfc::DFCalculation{QE}, threshold::Real; Epad=5.0)

Will set the energy window limits correctly according to the projections specified in the
structure of the job. The output of `projwfc` and the `threshold` will be used to determine
the minimum limit of the frozen energy window such that the interesting DOS of inside it exceeds
the threshold. `nscf` will be used to determine the DOS, and what the upper limit of the frozen window
needs to be to fit enough bands inside it, depending on the projections.
"""
function set_wanenergies!(job::DFJob, nscf::DFCalculation, projwfc::DFCalculation,
                          threshold::Real; Epad = 5.0)
    hasoutput_assert(projwfc)
    @assert isprojwfc(projwfc) "Please specify a valid projwfc calculation."
    @assert isnscf(nscf) "Please specify a valid nscf calculation."
    Emin = Emin_from_projwfc(job.structure, projwfc, threshold)
    return set_wanenergies!(job, nscf, Emin; Epad = Epad)
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
    nscf_calculation = getfirst(isnscf, calculations(job))
    projwfc_calculation = getfirst(isprojwfc, calculations(job))
    if projwfc_calculation === nothing || !hasoutput(projwfc_calculation)
        @info "No projwfc calculation found with valid output, using $min_window_determinator as Emin"
        return set_wanenergies!(job, nscf_calculation, min_window_determinator; kwargs...)
    else
        @info "Valid projwfc output found, using $min_window_determinator as the dos threshold."
        return set_wanenergies!(job, nscf_calculation, projwfc_calculation,
                                min_window_determinator; kwargs...)
    end
end

#--------------- Interacting with the Structure inside the DFJob ---------------#
#automatically sets the cell parameters for the entire job, implement others
Base.joinpath(job::DFJob, p) = joinpath(job.local_dir, p)

"""
    bandgap(job::DFJob, fermi=nothing)

Calculates the bandgap (possibly indirect) around the fermi level.
Uses the first found bands calculation, if there is none it uses the first found nscf calculation.
"""
function bandgap(job::DFJob, fermi = nothing)
    band_calcs = getfirst.([isbands, isnscf, isscf], (calculations(job),))
    if all(x -> x === nothing, band_calcs)
        error("No valid calculation found to calculate the bandgap.\nMake sure the job has either a valid bands or nscf calculation.")
    end
    if fermi === nothing
        fermi_calcs = getfirst.([isnscf, isscf], (calculations(job),))
        if all(x -> x === nothing, band_calcs)
            error("No valid calculation found to extract the fermi level.\nPlease supply the fermi level manually.")
        end
        fermi = maximum(readfermi.(filter(x -> x !== nothing, fermi_calcs)))
    end

    bands = readbands.(filter(x -> x !== nothing, band_calcs))
    return minimum(bandgap.(bands, fermi))
end

function readfermi(job::DFJob)
    ins = filter(x -> (isscf(x) || isnscf(x)) && hasoutfile(x), calculations(job))
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

function readbands(job::DFJob)
    calculation = getfirst(x -> isbands(x) && hasoutfile(x), calculations(job))
    if calculation === nothing
        calculation = getfirst(x -> isnscf(x) && hasoutfile(x), calculations(job))
        if calculation === nothing
            @warn "Job does not have a valid bands output."
            return nothing
        end
        @warn "No bands calculation found, return bands from nscf calculation."
        return readbands(calculation)
    end
    return readbands(calculation)
end

"""
    gencalc_wan(job::DFJob, min_window_determinator::Real, extra_wan_flags...; kwargs...)

Automates the generation of wannier calculations based on the `job`.
When a projwfc calculation is present in the `job`, `min_window_determinator` will be used to
determine the threshold value for including a band in the window based on the projections, otherwise
it will be used as the `Emin` value from which to start counting the number of bands needed for all
projections.
`extra_wan_flags` can be any extra flags for the Wannier90 calculation such as `write_hr` etc.
"""
function gencalc_wan(job::DFJob, min_window_determinator::Real, extra_wan_flags...;
                     kwargs...)
    nscf_calculation = getfirst(x -> isnscf(x), calculations(job))
    projwfc_calculation = getfirst(x -> isprojwfc(x), calculations(job))
    if projwfc_calculation === nothing || !hasoutput(projwfc_calculation)
        @info "No projwfc calculation found with valid output, using $min_window_determinator as Emin"
        return gencalc_wan(nscf_calculation, job.structure, min_window_determinator,
                           extra_wan_flags...; kwargs...)
    else
        @info "Valid projwfc output found, using $min_window_determinator as the dos threshold."
        return gencalc_wan(nscf_calculation, job.structure, projwfc_calculation,
                           min_window_determinator, extra_wan_flags...; kwargs...)
    end
end

#TODO: only for QE 
"Reads the pdos for a particular atom. Only works for QE."
function pdos(job::DFJob, atsym::Symbol, filter_word = "")
    projwfc = getfirst(isprojwfc, calculations(job))
    ats = atoms(job, atsym)
    @assert length(ats) > 0 "No atoms found with name $atsym."
    scf = getfirst(isscf, calculations(job))
    magnetic = any(ismagnetic, atoms(job)) || ismagnetic(scf)
    soc = issoc(scf)
    return pdos(projwfc, atsym, magnetic, soc, filter_word)
end

pdos(job::DFJob, atom::AbstractAtom, args...) = pdos(job, name(atom), args...)

function pdos(job::DFJob, atoms::Vector{AbstractAtom} = atoms(job), args...)
    t_energies, t_pdos = pdos(job, atoms[1], args...)
    for i in 2:length(atoms)
        t1, t2 = pdos(job, atoms[i], args...)
        t_pdos .+= t2
    end
    return (energies = t_energies, pdos = t_pdos)
end

"""
    last_submission(job::DFJob)

If a job was ever submitted, the last submission date is returned.
Otherwise 0 date is returned.
"""
function last_submission(job::DFJob)
    return get(job.metadata, :timestap, DateTime(0))
end

"""
    set_present!(job::DFJob, func)

Sets a function with the call signature `func(job)` which can be later called using the [`@present`](@ref) macro.
"""
function set_present!(job::DFJob, func)
    open(joinpath(job, ".present.jl"), "w") do f
        write(f, @code_string func(job))
    end
end

"""
    present(job)

Calls a present function if it was previously saved using [`set_present!`](@ref) or [`archive`](@ref). 
"""
macro present(job)
    return esc(quote
        if ispath(joinpath(job, ".present.jl"))
            t = include(joinpath(job, ".present.jl"))
            t(job)
        else
            @error "No presentation function defined.\n Please set it with `set_present!`."
        end
    end)
end

"""
    archive(job::DFJob, archive_directory::AbstractString, description::String=""; present = nothing, version=job.version)

Archives `job` by copying it's contents to `archive_directory` alongside a `results.jld2` file with all the parseable results as a Dict. `description` will be saved in a `description.txt` file in the `archive_directory`. A different job version can be copied using the `version` kwarg, and with the `present` kwarg a function can be specified that can be later called with the [`@present`](@ref) macro.
"""
function archive(job::DFJob, archive_directory::AbstractString, description::String=""; present = nothing, version=job.version)
    @assert !isarchived(job) "Job was already archived"
    final_dir = config_path(".archived", archive_directory)
    @assert !ispath(final_dir) "A archived job already exists in $archive_directory"
    mkpath(final_dir)

    out = outputdata(job)
    tj = deepcopy(job)
    switch_version!(tj, version)
    cp(tj, final_dir)
    set_localdir!(tj, final_dir)

    JLD2.save(joinpath(final_dir, "results.jld2"), "outputdata", out)
    
    present !== nothing && set_present!(tj, present)
    !isempty(description) && write(joinpath(final_dir, "description.txt"), description)
end
