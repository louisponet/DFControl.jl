const VERSION_DIR_NAME = ".versions"
const TEMP_CALC_DIR = "outputs"

@enum JobState BootFail Pending Running Completed Cancelled Deadline Failed NodeFail OutOfMemory Preempted Requeued Resizing Revoked Suspended Timeout Submitted Unknown 

"""
    Job(name::String, structure::Structure;
          calculations      ::Vector{Calculation} = Calculation[],
          dir               ::String = pwd(),
          header            ::Vector{String} = getdefault_jobheader(),
          metadata          ::Dict = Dict(),
          version           ::Int = last_job_version(dir),
          copy_temp_folders ::Bool = false, 
          server            ::String = getdefault_server(),
          environment::String ="")

A [`Job`](@ref) embodies a set of [`Calculations`](@ref Calculation) to be ran in directory `dir`, with the [`Structure`](@ref) as the subject.
## Keywords/further attributes
- `calculations`: calculations to calculations that will be run sequentially.
- `dir`: the directory where the calculations will be run.
- `header`: lines that will be pasted at the head of the job script, e.g. exports `export OMP_NUM_THREADS=1`, slurm settings`#SBATCH`, etc.
- `metadata`: various additional information, will be saved in `.metadata.jld2` in the `dir`.
- `version`: the current version of the job.
- `copy_temp_folders`: whether or not the temporary directory associated with intermediate calculation results should be copied when storing a job version. *CAUTION* These can be quite large.
- `server`: [`Server`](@ref) where to run the [`Job`](@ref).
- `environment`: [`Environment`](@ref) to be used for running the [`Job`](@ref).
 
    Job(job_name::String, structure::Structure, calculations::Vector{<:Calculation}, common_flags::Pair{Symbol, <:Any}...; kwargs...)

Creates a new job. The common flags will be attempted to be set in each of the `calculations`. The `kwargs...` are passed to the [`Job`](@ref) constructor. 

    Job(job_dir::String, job_script="job.tt"; version=nothing, kwargs...)

Loads the job in the `dir`.
If `job_dir` is not a valid job path, the previously saved jobs will be scanned for a job with a `dir` that
partly includes `job_dir`. If `version` is specified the corresponding job version will be returned if it exists. 
The `kwargs...` will be passed to the [`Job`](@ref) constructor.
"""
@with_kw_noshow mutable struct Job
    name::String = ""
    structure::Structure = Structure()
    calculations::Vector{Calculation} = Calculation[]
    dir::String = pwd()
    header::Vector{String} = String[]
    metadata::Dict{Symbol,Any} = Dict{Symbol,Any}()
    version::Int = -1
    copy_temp_folders::Bool = false
    server::String = "localhost"
    environment::String = ""
    function Job(name, structure, calculations, dir, header, metadata, version,
                 copy_temp_folders, server, environment)
        if dir[end] == '/'
            dir = dir[1:end-1]
        end
        out = new(name, structure, calculations, dir, header, metadata, version,
                  copy_temp_folders, server, environment)
        return out
    end
end

#TODO implement abinit
function Job(job_name::String, structure::Structure, calculations::Vector{<:Calculation},
             common_flags::Pair{Symbol,<:Any}...; kwargs...)
    out = Job(; name = job_name, structure = structure, calculations = calculations,
              kwargs...)
    for (f, v) in common_flags
        out[f] = v
    end
    return out
end

StructTypes.StructType(::Type{Job}) = StructTypes.Mutable()

#-------------------BEGINNING GENERAL SECTION-------------#
scriptpath(job::Job) = joinpath(job, "job.tt")
submission_time(job::Job)  = mtime(scriptpath(job))

runslocal(job::Job) = job.server == "localhost"
isarchived(job::Job) = occursin(".archived", job.dir)

Base.abspath(job::Job) = abspath(Server(job.server), job.dir)
    
Base.ispath(job::Job, p...) =
    runslocal(job) ? joinpath(job, p...) : ispath(DFC.Servers.maybe_start_server(DFC.Server(job.server)), joinpath(job, p...))

    
"""
    joinpath(job::Job, args...)

`joinpath(job.dir, args...)`.
"""
Base.joinpath(job::Job, args...) = joinpath(abspath(job), args...)
Base.readdir(job::Job) = runslocal(job) ? readdir(abspath(job)) : readdir(Server(job.server), job.dir)
    
function Base.pop!(job::Job, name::String)
    i = findfirst(x -> x.name == name, job.calculations)
    if i === nothing
        error("No calculation with name $name found.")
    end
    out = job.calculations[i]
    deleteat!(job.calculations, i)
    return out
end

function Utils.searchdir(job::Job, str::AbstractString)
    if runslocal(job)
        return searchdir(abspath(job), str)
    else
        return searchdir(Server(job.server), abspath(job), str)
    end
end

function Base.setindex!(job::Job, value, key::Symbol)
    for c in job.calculations
        c[key] = value
    end
end

function Base.getindex(job::Job, flg::Symbol)
    outdict = Dict()
    for i in job.calculations
        tfl = get(i, flg, nothing)
        if tfl !== nothing
            outdict[i.name] = tfl
        end
    end
    return outdict
end

function Base.pop!(job::Job, args...)
    out = Dict{String, Any}()
    for c in job.calculations
        out[c.name] = pop!(c, args...)
    end
    return out
end

"""
    insert!(job::Job, i::Int, calculation::Calculation) = insert!(job.calculations, i, calculation)
"""
function Base.insert!(job::Job, index::Int, calculation::Calculation)
    return insert!(job.calculations, index, calculation)
end

"""
    push!(job::Job, calculation::Calculation) = push!(job.calculations, calculation)
"""
Base.push!(job::Job, calculation::Calculation) = push!(job.calculations, calculation)

"""
    pop!(job::Job) = pop!(job.calculations)
"""
Base.pop!(job::Job) = pop!(job.calculations)

"""
    append!(job::Job, args...) = append!(job.calculations, args...)
"""
Base.append!(job::Job, args...) = append!(job.calculations, args...)
Base.length(job::Job) = length(job.calculations)
Base.lastindex(job::Job) = length(job)
Base.getindex(job::Job, i::Integer) = job.calculations[i]

"""
    getindex(job::Job, name::String)
    
Returns the `Calculation` with the specified `name`.

    getindex(job::Job, i::Integer)
    
Returns the i'th `Calculation` in the job.
"""
function Base.getindex(job::Job, id::String)
    tmp = getfirst(x -> x.name == id, job.calculations)
    if tmp != nothing
        return tmp
    else
        error("Calculation $id not found.")
    end
end

Base.getindex(job::Job, el::Element) = job.structure[el]

"""
    main_job_dir(dir::AbstractString)
    main_job_dir(job::Job)

Returns the main directory of the job, also when the job's version is not the one
in the main directory.
"""
function main_job_dir(dir::AbstractString)
    d = split(dir, Jobs.VERSION_DIR_NAME)[1]
    return d[end] == '/' ? d[1:end-1] : d
end
main_job_dir(job::Job) = isabspath(job.dir) ? main_job_dir(job.dir) : main_job_dir(joinpath(Server(job), job.dir))

"""
    gencalc_wan(job::Job, min_window_determinator::Real, extra_wan_flags...; kwargs...)

Automates the generation of wannier calculations based on the `job`.
When a projwfc calculation is present in the `job`, `min_window_determinator` will be used to
determine the threshold value for including a band in the window based on the projections, otherwise
it will be used as the `Emin` value from which to start counting the number of bands needed for all
projections.
`extra_wan_flags` can be any extra flags for the Wannier90 calculation such as `write_hr` etc.
"""
function Calculations.gencalc_wan(job::Job, min_window_determinator::Real,
                                  extra_wan_flags...; kwargs...)
    nscf_calc = getfirst(x -> Calculations.isnscf(x), job.calculations)
    nscf_calc === nothing && "Please first run an nscf calculation."
    if get(nscf_calc, :nosym, false) != true
        @info "'nosym' flag was not set in the nscf calculation.
                If this was not intended please set it and rerun the nscf calculation.
                This generally gives errors because of omitted kpoints, needed for pw2wannier90.x"
    end
    projwfc_calc = getfirst(x -> Calculations.isprojwfc(x), job.calculations)
    outdat = outputdata(job)
    if projwfc_calc === nothing || !haskey(outdat, projwfc_calc.name)
        @info "No projwfc calculation found with valid output, using $min_window_determinator as Emin"
        Emin = min_window_determinator
    else
        @info "Valid projwfc output found, using $min_window_determinator as the dos threshold."
        Emin = Calculations.Emin_from_projwfc(job.structure, outdat[projwfc_calc.name][:states], outdat[projwfc_calc.name][:bands], min_window_determinator)
    end
    return Calculations.gencalc_wan(nscf_calc, job.structure, outdat[nscf_calc.name][:bands], Emin, extra_wan_flags...; kwargs...)
end

"""
    set_flow!(job::Job, should_runs::Pair{String, Bool}...)

Sets whether or not calculations should be scheduled to run.
The `name` of each calculation in the job will be checked against the string in each pair of `should_runs`, and the
`calculation.run` will be set accordingly.

Example:
```julia
set_flow!(job, "" => false, "scf" => true)
```
would un-schedule all calculations in the job, and schedule the "scf" and "nscf" calculations to run.
"""
function set_flow!(job::Job, should_runs...)
    for (name, run) in should_runs
        for calculation in filter(x -> occursin(name, x.name), job.calculations)
            calculation.run = run
        end
    end
    return job
end

"""
    set_headerword!(job::Job, old_new::Pair{String, String})

Replaces the specified word in the header with the new word.
"""
function set_headerword!(job::Job, old_new::Pair{String,String}; print = true)
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

#TODO
"""
    set_wanenergies!(job::Job, nscf::Calculation{QE}, projwfc::Calculation{QE}, threshold::Real; Epad=5.0)

Will set the energy window limits correctly according to the projections specified in the
structure of the job. The output of `projwfc` and the `threshold` will be used to determine
the minimum limit of the frozen energy window such that the interesting DOS of inside it exceeds
the threshold. `nscf` will be used to determine the DOS, and what the upper limit of the frozen window
needs to be to fit enough bands inside it, depending on the projections.
"""
function set_wanenergies!(job::Job, nscf::Calculation, projwfc::Calculation,
                          threshold::Real; Epad = 5.0)
    hasoutput_assert(projwfc)
    @assert Calculations.isprojwfc(projwfc) "Please specify a valid projwfc calculation."
    @assert Calculations.isnscf(nscf) "Please specify a valid nscf calculation."
    Emin = Emin_from_projwfc(job.structure, projwfc, threshold)
    return set_wanenergies!(job, nscf, Emin; Epad = Epad)
end

"""
    set_wanenergies!(job::Job, min_window_determinator::Real; kwargs...)

Sets the energy windows of wannier calculations based on the `job`.
When a projwfc calculation is present in the `job`, `min_window_determinator` will be used to
determine the threshold value for including a band in the window based on the projections, otherwise
it will be used as the `Emin` value from which to start counting the number of bands needed for all
projections.
"""
function set_wanenergies!(job::Job, min_window_determinator::Real; kwargs...)
    nscf_calc = getfirst(Calculations.isnscf, job.calculations)
    projwfc_calc = getfirst(Calculations.isprojwfc, job.calculations)
    if projwfc_calc === nothing || !hasoutput(projwfc_calculation)
        @info "No projwfc calculation found with valid output, using $min_window_determinator as Emin"
        return set_wanenergies!(job, nscf_calc, min_window_determinator; kwargs...)
    else
        @info "Valid projwfc output found, using $min_window_determinator as the dos threshold."
        return set_wanenergies!(job, nscf_calc, projwfc_calculation,
                                min_window_determinator; kwargs...)
    end
end

"""
    bandgap(job::Job, fermi=nothing)

Calculates the bandgap (possibly indirect) around the fermi level.
Uses the first found bands calculation, if there is none it uses the first found nscf calculation.
"""
function bandgap(job::Job, fermi = nothing)
    band_calcs = filter(!isnothing, getfirst.([Calculations.isbands, Calculations.isnscf, Calculations.isscf], (job.calculations,)))
    if isempty(band_calcs)
        error("No valid calculation found to calculate the bandgap.\nMake sure the job has either a valid bands or nscf calculation.")
    end

    outdat = outputdata(job)
    if fermi === nothing
        fermi_calcs = filter(!isnothing, getfirst.([Calculations.isnscf, Calculations.isscf], (job.calculations,)))
        if isempty(fermi_calcs)
            error("No valid calculation found to extract the fermi level.\nPlease supply the fermi level manually.")
        end
        fermi = maximum(x->get(get(outdat, x.name, Dict()), :fermi, -Inf), fermi_calcs)
    end

    bands = map(x->outdata[x.name][:bands], filter(x -> haskey(outdat, x.name), band_calcs))
    return minimum(bandgap.(bands, fermi))
end

function readfermi(job::Job, outdat)
    ins = filter(x -> (Calculations.isvcrelax(x) || Calculations.isscf(x) || Calculations.isnscf(x)), job.calculations)
    @assert isempty(ins) !== nothing "Job does not have a valid scf or nscf output."
    for i in ins
        if haskey(outdat, i.name)
            o = outdat[i.name]
            if haskey(o, :fermi)
                return o[:fermi]
            end
        end
    end
    @warn "No output files with fermi level found."
    return 0.0
end

function readbands(job::Job, outdat)
    calc = getfirst(x -> Calculations.isbands(x), job.calculations)
    outdat = outputdata(job)
    if calc === nothing || !haskey(outdat, calc.name)
        calc = getfirst(x -> Calculations.isnscf(x), job.calculations)
        if calc === nothing || !haskey(outdat, calc.name)
            @warn "Job does not have a valid bands output."
            return nothing
        end
        @warn "No bands calculation found, return bands from nscf calculation."
        return outdat[calc.name][:bands]
    end
    return outdat[calc.name][:bands]
end

#TODO: only for QE 
"""
    set_present!(job::Job, func::Function)
    set_present!(job::Job, func::String)
    set_present!(job::Job, func::Expr)

Sets a function with the call signature `func(job)` which can be later called using the [`@present`](@ref) macro.
"""
function set_present!(job::Job, func::Function)
    try
        str = loaded_modules_string() * @code_string func(job)
        set_present!(job, str)
    catch
        error("Could not generate the source string for the supplied function.\nIf you are running in a Jupyter notebook, please supply the source code as a string, or expression to set_present!.\nDon't forget to include the right using statements.")
    end
end
function set_present!(job::Job, func::AbstractString)
    open(joinpath(job, ".present.jl"), "w") do f
        return write(f, func)
    end
end
function set_present!(job::Job, func::Expr)
    funcstr = string(func)
    if funcstr[1:5] == "begin"
        funcstr = funcstr[findfirst(isequal('\n'), funcstr)+1:findlast(isequal('\n'), funcstr)-1]
    end
    return set_present!(job, funcstr)
end

"""
    present(job)

Calls a present function if it was previously saved using [`set_present!`](@ref) or [`archive`](@ref). 
"""
macro present(job)
    return esc(quote
                   if ispath(joinpath($job, ".present.jl"))
                       t = include(joinpath($job, ".present.jl"))
                       DFControl.with_logger(DFControl.MinLevelLogger(DFControl.current_logger(),
                                                                      DFControl.Logging.Error)) do
                           return t($job)
                       end
                   else
                       @error "No presentation function defined.\n Please set it with `set_present!`."
                   end
               end)
end

"""
    archive(job::Job, archive_directory::AbstractString, description::String=""; present = nothing, version=job.version)

Archives `job` by copying it's contents to `archive_directory` alongside a `results.jld2` file with all the parseable results as a Dict. `description` will be saved in a `description.txt` file in the `archive_directory`. A different job version can be copied using the `version` kwarg, and with the `present` kwarg a function can be specified that can be later called with the [`@present`](@ref) macro.
"""
function archive(job::Job, archive_directory::AbstractString, description::String = "";
                 present = nothing, version = job.version)
    @assert !isarchived(job) "Job was already archived"
    final_dir = config_path(".archived", archive_directory)
    @assert !ispath(final_dir) "A archived job already exists in $archive_directory"
    mkpath(final_dir)

    present !== nothing && set_present!(job, present)
    out = outputdata(job)
    tj = deepcopy(job)
    switch_version!(tj, version)
    cp(tj, final_dir)
    tj.dir = final_dir

    JLD2.jldsave(joinpath(final_dir, "results.jld2"); outputdata=out)

    !isempty(description) && write(joinpath(final_dir, "description.txt"), description)
    push!(JOB_REGISTRY.archived, tj.dir)
    @info "Archived job at $(tj.dir). If you're done with this one, it is safe to delete the directory at $(job.dir)."
    write_job_registry()
    return nothing
end

for (f, strs) in zip((:cp, :mv), (("copy", "Copies"), ("move", "Moves")))
    @eval begin
        """
            $($f)(job::Job, dest::AbstractString; all=false, temp=false, kwargs...)

        $($(strs[2])) the contents of `job.dir` to `dest`. If `all=true`, it will also $($(strs[1])) the
        `.version` directory with all previous versions. If `temp=true` it will override
        `job.copy_temp_folders` and $($(strs[1])) also the temporary calculation directories.
        The `kwargs...` are passed to `Base.$($f)`.
        """
        function Base.$f(job::Job, dest::AbstractString; all = false, temp = false,
                         kwargs...)
            if !ispath(dest)
                mkpath(dest)
            end
            for file in readdir(job)
                if !all
                    if file == VERSION_DIR_NAME
                        continue
                    elseif file == TEMP_CALC_DIR && !(temp || job.copy_temp_folders)
                        continue
                    end
                end
                if joinpath(job, file) == abspath(dest)
                    continue
                end
                $f(joinpath(job, file), joinpath(dest, file); kwargs...)
            end
        end
    end
end

function pseudo_cutoffs(str::String)
    lines = split(str, "\n")
    id = findfirst(x->occursin("Suggested minimum cutoff for wavefunctions:", x), lines)
    if id !== nothing
        return parse(Float64, split(lines[id])[end-1]), parse(Float64, split(lines[id+1])[end-1])
    else
        return 0.0, 0.0
    end
end

function sanitize_cutoffs!(job::Job)
    ψcut, ρcut = 0.0, 0.0
    # the assumption is that the most important cutoff calculation is the scf/vcrelax that is ran first 
    ψ_cut_calc = getfirst(x -> haskey(x, Calculations.ψ_cutoff_flag(x)), job.calculations)
    if ψ_cut_calc !== nothing
        ψcut = ψ_cut_calc[Calculations.ψ_cutoff_flag(ψ_cut_calc)]
    else
        pseudos = unique(map(x -> x.pseudo, job.structure.atoms))
        all_cuts = map(x -> pseudo_cutoffs(x), pseudos) 
        ψcut = maximum(x -> x[1], all_cuts)
        ρcut = maximum(x -> x[2], all_cuts)
        @assert ψcut != 0.0 "No energy cutoff was specified in any calculation, and the calculated cutoff from the pseudopotentials was 0.0.\nPlease manually set one."
        @info "No energy cutoff was specified in the scf calculation.\nCalculated ψcut=$ψcut."
    end
    ρ_cut_calc = getfirst(x -> Calculations.hasflag(x, Calculations.ρ_cutoff_flag(x)),
                          job.calculations)
    if ρ_cut_calc !== nothing
        ρcut = ρ_cut_calc[Calculations.ρ_cutoff_flag(ρ_cut_calc)]
    end
    for i in job.calculations
        ψflag = Calculations.ψ_cutoff_flag(i)
        ψflag !== nothing &&
            !haskey(i, ψflag) &&
            Calculations.set_flags!(i, ψflag => ψcut; print = false)
        ρflag = Calculations.ρ_cutoff_flag(i)
        ρflag !== nothing && !haskey(i, ρflag) && ρcut != 0.0 &&
            Calculations.set_flags!(i, ρflag => ρcut; print = false)
    end
end

