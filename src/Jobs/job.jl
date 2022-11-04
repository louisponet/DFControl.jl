const VERSION_DIR_NAME = ".versions"
const TEMP_CALC_DIR = "outputs"

@enum JobState BootFail Pending Running Completed Configuring Completing Cancelled Deadline Failed NodeFail OutOfMemory Preempted Requeued Resizing Revoked Suspended Timeout Submitted Unknown PostProcessing Saved

"""
    Job(name::String, structure::Structure;
          calculations      ::Vector{Calculation} = Calculation[],
          dir               ::String = pwd(),
          version           ::Int = last_job_version(dir),
          copy_temp_folders ::Bool = false, 
          server            ::String = getdefault_server(),
          environment::String ="")

A [`Job`](@ref) embodies a set of [`Calculations`](@ref Calculation) to be ran in directory `dir`, with the [`Structure`](@ref) as the subject.
## Keywords/further attributes
- `calculations`: calculations to calculations that will be run sequentially.
- `dir`: the directory where the calculations will be run.
- `version`: the current version of the job.
- `copy_temp_folders`: whether or not the temporary directory associated with intermediate calculation results should be copied when storing a job version. *CAUTION* These can be quite large.
- `server`: [`Server`](@ref) where to run the [`Job`](@ref).
- `environment`: [`Environment`](@ref) to be used for running the [`Job`](@ref).
 
    Job(job_name::String, structure::Structure, calculations::Vector{<:Calculation}, common_flags::Pair{Symbol, <:Any}...; kwargs...)

Creates a new job. The common flags will be attempted to be set in each of the `calculations`. The `kwargs...` are passed to the [`Job`](@ref) constructor. 

    Job(job_dir::String, job_script="job.sh"; version=nothing, kwargs...)

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
    version::Int = 0
    copy_temp_folders::Bool = false
    server::String = gethostname()
    environment::String = ""
    function Job(name, structure, calculations, dir, version,
                 copy_temp_folders, server, environment)
        if !isempty(dir)
            if dir[end] == '/'
                dir = dir[1:end-1]
            end
        end
        out = new(name, structure, calculations, dir, version,
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
Job(dir::AbstractString; kwargs...) = Job(;dir=dir, kwargs...)

StructTypes.StructType(::Type{Job}) = StructTypes.Mutable()

#-------------------BEGINNING GENERAL SECTION-------------#
scriptpath(job::Job) = joinpath(job, "job.sh")
submission_time(job::Job)  = mtime(scriptpath(job))

runslocal(job::Job) = job.server == "localhost"
isarchived(job::Job) = occursin("archived", job.dir)

"""
    abspath(job::Job, args...)

If the job is local this is `abspath(job.dir)`, otherwise it will resolve the abspath using the [`Server`](@ref) rootdir.
"""
Base.abspath(job::Job) = abspath(Server(job.server), job.dir)
    
Base.ispath(job::Job, p...) =
    runslocal(job) ? joinpath(job, p...) : ispath(Server(job.server), joinpath(job, p...))

    
"""
    joinpath(job::Job, args...)

If the job is local this is `joinpath(job.dir, args...)`, otherwise it will resolve the path using the [`Server`](@ref) rootdir.
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

Utils.searchdir(job::Job, str::AbstractString) = searchdir(Server(job.server), abspath(job), str)

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
    if tmp !== nothing
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
    if projwfc_calc === nothing
        @info "No projwfc calculation found with valid output, using $min_window_determinator as Emin"
        return set_wanenergies!(job, nscf_calc, min_window_determinator; kwargs...)
    else
        @info "Valid projwfc output found, using $min_window_determinator as the dos threshold."
        return set_wanenergies!(job, nscf_calc, projwfc_calc,
                                min_window_determinator; kwargs...)
    end
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
        function Base.$f(job::Job, dest::AbstractString; 
                         kwargs...)
            server = Server(job.server)
            if !ispath(server, dest)
                mkpath(server, dest)
            end
            for file in readdir(server, job.dir)
                p = joinpath(job, file)
                if file == VERSION_DIR_NAME
                    continue
                elseif file == TEMP_CALC_DIR && !(temp || job.copy_temp_folders)
                    continue
                end
                if joinpath(job, file) == abspath(dest)
                    continue
                end
                $f(server, p, joinpath(dest, file); kwargs...)
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
    ψ_cut_calc = getfirst(x -> Calculations.ψ_cutoff_flag(x) !== nothing && haskey(x, Calculations.ψ_cutoff_flag(x)), job.calculations)
    if ψ_cut_calc !== nothing
        ψcut = ψ_cut_calc[Calculations.ψ_cutoff_flag(ψ_cut_calc)]
    else
        pseudos = unique(map(x -> x.pseudo, job.structure.atoms))
        pseudo_strings = map(pseudos) do ps
            if !isempty(ps.pseudo)
                return ps.pseudo
            elseif !isempty(ps.server)
                s = Server(ps.server)
                if isalive(s) && ispath(s, ps.path)
                    return read(s, ps.path, String)
                else
                    return ""
                end
            end
        end
        if !isempty(pseudo_strings) && !all(isequal(nothing), pseudo_strings)
            all_cuts = map(x -> pseudo_cutoffs(x), pseudo_strings) 
            ψcut = maximum(x -> x[1], all_cuts)
            ρcut = maximum(x -> x[2], all_cuts)
            @assert ψcut != 0.0 "No energy cutoff was specified in any calculation, and the calculated cutoff from the pseudopotentials was 0.0.\nPlease manually set one."
            @info "No energy cutoff was specified in the scf calculation.\nCalculated ψcut=$ψcut."
        else
            return @warn "No cutoffs found or set."
        end
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

function RemoteHPC.pull(j::Job, f, t)
    @assert ispath(j, f) "File $f not found in jobdir."
    RemoteHPC.pull(Server(j.server), joinpath(j, f), t)
end

function timestamp(jobdir::AbstractString)
    scriptpath = joinpath(jobdir, "job.sh")
    if ispath(scriptpath)
        return unix2datetime(mtime(scriptpath))
    else
        return DateTime(0)
    end
end


