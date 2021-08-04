const VERSION_DIR_NAME = ".versions"
const TEMP_CALC_DIR = "outputs"

"""
    Job(name::String, structure::Structure;
          calculations      ::Vector{Calculation} = Calculation[],
          dir               ::String = pwd(),
          header            ::Vector{String} = getdefault_jobheader(),
          metadata          ::Dict = Dict(),
          version           ::Int = last_job_version(dir),
          copy_temp_folders ::Bool = false, 
          server            ::String = getdefault_server())

A [`Job`](@ref) embodies a set of [`Calculations`](@ref Calculation) to be ran in directory `dir`, with the [`Structure`](@ref) as the subject.
## Keywords/further attributes
- `calculations`: calculations to calculations that will be run sequentially.
- `dir`: the directory where the calculations will be run.
- `header`: lines that will be pasted at the head of the job script, e.g. exports `export OMP_NUM_THREADS=1`, slurm settings`#SBATCH`, etc.
- `metadata`: various additional information, will be saved in `.metadata.jld2` in the `dir`.
- `version`: the current version of the job.
- `copy_temp_folders`: whether or not the temporary directory associated with intermediate calculation results should be copied when storing a job version. *CAUTION* These can be quite large.
- `server`: [`Server`](@ref) where to run the [`Job`](@ref).
 
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
    function Job(name, structure, calculations, dir, header, metadata, version,
                   copy_temp_folders, server)
        if dir[end] == '/'
            dir = dir[1:end-1]
        end
        if !isabspath(dir)
            dir = abspath(dir)
        end
        if isempty(structure.name)
            structure.name = split(name, "_")[1]
        end
        if isempty(metadata)
            mpath = joinpath(dir, ".metadata.jld2")
            if ispath(mpath)
                stored_data = JLD2.load(mpath)
                metadata = haskey(stored_data, "metadata") ? stored_data["metadata"] :
                           metadata
                version = haskey(stored_data, "version") ? stored_data["version"] : version
            end
        end
        out = new(name, structure, calculations, dir, header, metadata, version,
                  copy_temp_folders, server)
        return out
    end
end

#TODO implement abinit
function Job(job_name::String, structure::Structure,
               calculations::Vector{<:Calculation}, common_flags::Pair{Symbol,<:Any}...;
               kwargs...)
    out = Job(; name = job_name, structure = structure, calculations = calculations,
                kwargs...)
    for (f, v) in common_flags
        out[f] = v
    end
    return out
end

StructTypes.StructType(::Type{Job}) = StructTypes.Mutable()


#-------------------BEGINNING GENERAL SECTION-------------#
scriptpath(job::Job) = joinpath(job.dir, "job.tt")
starttime(job::Job)  = mtime(scriptpath(job))

runslocal(job::Job)    = job.server == "localhost"
isarchived(job::Job) = occursin(".archived", job.dir)

"""
    joinpath(job::Job, args...)

`joinpath(job.dir, args...)`.
"""
Base.joinpath(job::Job, args...) = joinpath(job.dir, args...)

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
    return joinpath.((job,), searchdir(job.dir, str))
end

function Base.setindex!(job::Job, value, key::Symbol)
    for c in job.calculations
        c[key] = value
    end
end

function Base.getindex(job::Job, flg::Symbol)
    outdict = Dict()
    for i in job.calculations
        tfl = i[flg]
        if tfl !== nothing
            outdict[i.name]= tfl
        end
    end
    return outdict
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
    gencalc_wan(job::Job, min_window_determinator::Real, extra_wan_flags...; kwargs...)

Automates the generation of wannier calculations based on the `job`.
When a projwfc calculation is present in the `job`, `min_window_determinator` will be used to
determine the threshold value for including a band in the window based on the projections, otherwise
it will be used as the `Emin` value from which to start counting the number of bands needed for all
projections.
`extra_wan_flags` can be any extra flags for the Wannier90 calculation such as `write_hr` etc.
"""
function Calculations.gencalc_wan(job::Job, min_window_determinator::Real, extra_wan_flags...;
                     kwargs...)
    nscf_calculation = getfirst(x -> isnscf(x), job.calculations)
    projwfc_calculation = getfirst(x -> isprojwfc(x), job.calculations)
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

"""
    main_job_dir(dir::AbstractString)
    main_job_dir(job::Job)

Returns the main directory of the job, also when the job's version is not the one
in the main directory.
"""
main_job_dir(dir::AbstractString) = split(dir, Jobs.VERSION_DIR_NAME)[1]
main_job_dir(job::Job) = main_job_dir(job.dir)
