const VERSION_DIR_NAME = ".versions"
const TEMP_CALC_DIR = "outputs"

name(job) = job.name
#-------------------BEGINNING GENERAL SECTION-------------#
scriptpath(job::DFJob) = joinpath(job.dir, "job.tt")
starttime(job::DFJob)  = mtime(scriptpath(job))

runslocal(job::DFJob)    = job.server == "localhost"
structure(job::DFJob)    = job.structure
isQEjob(job::DFJob)      = any(x -> package(x) == QE, job.calculations)
iswannierjob(job::DFJob) = any(x -> package(x) == Wannier90, job.calculations) && any(x -> isnscf(x), job.calculations)
getnscfcalc(job::DFJob)  = getfirst(x -> isnscf(x), job.calculations)

cell(job::DFJob)         = cell(structure(job))
calculations(job::DFJob) = job.calculations
isarchived(job::DFJob) = occursin(".archived", job.dir)


"""
    joinpath(job::DFJob, args...)

`joinpath(job.dir, args...)`.
"""
Base.joinpath(job::DFJob, args...) = joinpath(job.dir, args...)

function Base.pop!(job::DFJob, name::String)
    i = findfirst(x -> x.name == name, job.calculations)
    if i === nothing
        error("No calculation with name $name found.")
    end
    out = job.calculations[i]
    deleteat!(job.calculations, i)
    return out
end

function Utils.searchdir(job::DFJob, str::AbstractString)
    return joinpath.((job,), searchdir(job.dir, str))
end

function Base.getindex(job::DFJob, flg::Symbol)
    outdict = Dict()
    for i in calculations(job)
        tfl = flag(i, flg)
        if tfl != nothing
            outdict[name(i)] = tfl
        end
    end
    return outdict
end

function set_flags!(job::DFJob, calculations::Vector{<:DFCalculation}, flags...;
                    print = true)
    found_keys = Symbol[]

    for calc in calculations
        t_, = set_flags!(calc, flags...; print = print)
        push!(found_keys, t_...)
    end
    nfound = setdiff([k for (k, v) in flags], found_keys)
    if print && length(nfound) > 0
        f = length(nfound) == 1 ? "flag" : "flags"
        dfprintln("$f '$(join(":" .* String.(setdiff(flagkeys, found_keys)),", "))' were not found in the allowed calculation variables of the specified calculations!")
    end
    return job
end

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
        error("Calculation $id not found.")
    end
end

Base.getindex(job::DFJob, el::Element) = job.structure[el]

