using ..Utils: string2cmd

abstract type Scheduler end
@with_kw struct Slurm <: Scheduler
    type::String = "slurm"
end
@with_kw struct Bash <: Scheduler
    type::String="bash"
end
@with_kw struct HQ <: Scheduler
    type::String="hq"
    server_command::String="hq"
end
StructTypes.StructType(::Type{Cmd}) = StructTypes.Struct()
StructTypes.StructType(::Type{Scheduler}) = StructTypes.AbstractType()
StructTypes.subtypes(::Type{Scheduler}) = (bash = Bash, slurm = Slurm, hq=HQ)

StructTypes.StructType(::Type{Bash})  = StructTypes.Struct()
StructTypes.StructType(::Type{Slurm}) = StructTypes.Struct()
StructTypes.StructType(::Type{HQ})    = StructTypes.Struct()
StructTypes.subtypekey(::Type{Scheduler}) = :type

submit_cmd(s::S) where {S<:Scheduler} = error("No submit_cmd method defined for $S.")
submit_cmd(s::Slurm) = `sbatch`
submit_cmd(s::Bash)  = `bash`
submit_cmd(s::HQ)  = `hq`

function is_reachable(server_command::String)
    t = run(string2cmd("which $server_command"), wait=false)
    while !process_exited(t)
        sleep(0.005)
    end
    return t.exitcode == 0
end
is_reachable(s::HQ) = is_reachable(s.server_command)
is_reachable(::Slurm) = is_reachable("sbatch")
is_reachable(::Bash) = true


scheduler_directive_prefix(::Slurm) = "#SBATCH"
scheduler_directive_prefix(::Bash) = "#"
scheduler_directive_prefix(::HQ) = "#HQ"

scheduler_name_flag(::Slurm) = "job-name"
scheduler_name_flag(::Bash) = "job-name"
scheduler_name_flag(::HQ) = "name"

# These will be filled in by definitions in Service
submit(s::S, jobdir::String) where {S<:Scheduler} = error("No submit method defined for $S.")
abort(s::S, id::Int) where {S<:Scheduler} = error("No abort method defined for $S.")
jobstate(s::S, id) where {S<:Scheduler} = error("No jobstate method defined for $S.")
jobid(s::S, dir::AbstractString) where {S<:Scheduler} = error("No jobid method defined for $S.")
