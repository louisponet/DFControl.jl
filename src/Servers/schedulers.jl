abstract type Scheduler end
@with_kw struct Slurm <: Scheduler
    type::String = "slurm"
end
@with_kw struct Bash <: Scheduler
    type::String="bash"
end
@with_kw struct HQ <: Scheduler
    type::String="hq"
end
StructTypes.StructType(::Type{Scheduler}) = StructTypes.AbstractType()
StructTypes.subtypes(::Type{Scheduler}) = (bash = Bash, slurm = Slurm, hq=HQ)

StructTypes.StructType(::Type{Bash})  = StructTypes.Struct()
StructTypes.StructType(::Type{Slurm}) = StructTypes.Struct()
StructTypes.StructType(::Type{HQ})    = StructTypes.Struct()
StructTypes.subtypekey(::Type{Scheduler}) = :type

submit_cmd(s::S) where {S<:Scheduler} = error("No submit_cmd method defined for $S.")
submit_cmd(s::Slurm) = `sbatch`
submit_cmd(s::Bash)  = `bash`
submit_cmd(s::Bash)  = `hq`

# These will be filled in by definitions in Service
submit(s::S, jobdir::String) where {S<:Scheduler} = error("No submit method defined for $S.")
abort(s::S, id::Int) where {S<:Scheduler} = error("No abort method defined for $S.")
jobstate(s::S, id) where {S<:Scheduler} = error("No jobstate method defined for $S.")
jobid(s::S, dir::AbstractString) where {S<:Scheduler} = error("No jobid method defined for $S.")
