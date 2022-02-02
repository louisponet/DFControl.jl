"""
    Environment(MPI_command::String, scheduler_flags::Vector{String}, exports::Vector{String})

Environment to run a [`Job`](@ref) in. When running on a server with a scheduler `scheduler_flags` holds what these should be. e.g. `#SBATCH -N 2`. `MPI_command` and `MPI_processes` will be used to prefix executables that should be ran in parallel, e.g. if `MPI_command = "mpirun -np 4 --bind-to core"` a parallel executable will be translated into a script line as `mpirun -np 4 --bind-to core exec`. 
"""
@with_kw mutable struct Environment <: Storable
    name::String = ""
    MPI_command::String = ""
    scheduler_flags::Vector{String} = String[]
    exports::Vector{String} = String[]
end
Environment(name::AbstractString; kwargs...) = Environment(;name=name, kwargs...)

StructTypes.StructType(::Type{Environment}) = StructTypes.Struct()

Database.storage_directory(e::Environment) = "environments"

function Base.:(==)(e1::Environment, e2::Environment)
    if e1.MPI_command != e2.MPI_command || e1.name != e2.name
        return false
    else
        for fn in (:scheduler_flags, :exports)
            for f1 in getfield(e1, fn)
                if !(f1 ∈ getfield(e2, fn))
                    return false
                end
            end
        end
    end
    return true
end

"Reads the script file and tries to understand what the runtime environment is."
function environment_from_jobscript(scriptpath::String)
    lines = filter(!isempty, strip.(readlines(scriptpath)))
    id = findfirst(l -> any(x -> occursin(x, l), Calculations.RUN_EXECS), lines)
    if id !== nothing
        sline = split(replace(lines[id], "#" => ""))

        MPI_command = join(sline[1:findnext(x -> splitdir(x)[end] ∈ Calculations.allexecs() , sline, 2)-1], " ")
    else
        MPI_command = ""
    end
    #TODO only works with slurm!!!
    scheduler_flags = map(y -> strip(split(y, "SBATCH")[end]), filter(x -> occursin("#SBATCH", x), lines))
    deleteat!(scheduler_flags, findall(x -> occursin("-J", x) || occursin("--job-name", x), scheduler_flags))
    exports = map(y -> strip(split(y, "export")[end]), filter(x -> occursin("export", x), lines))
    
    t = Environment("", MPI_command, scheduler_flags, exports)
    n = Database.name(t)
    t.name = n === nothing ? t.name : n
    return t
end

function Base.write(f::IO, env::Environment)
    for flag in env.scheduler_flags
        write(f, "#SBATCH $flag\n")
    end
    for flag in env.exports
        write(f, "export $flag\n")
    end
end

function Base.rm(e::Environment)
    fullpath = joinpath(ENVIRONMENTS_DIR, e.name * ".json")
    if ispath(fullpath)
        rm(fullpath)
    end
end
