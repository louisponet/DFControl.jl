#here all the different input structures for the different calculations go
#TODO should we also create a config file for each job with stuff like server etc? and other config things,
#      which if not supplied could contain the default stuff?
"""
Represents a full DFT job with multiple input files and calculations.
"""
@with_kw_noshow mutable struct DFJob
    name         ::String
    structure    ::AbstractStructure
    inputs       ::Vector{DFInput} = DFInput[]
    local_dir    ::String = pwd()
    server       ::String = getdefault_server()
    server_dir   ::String = ""
    header       ::Vector{String} = getdefault_jobheader()
    metadata     ::Dict = Dict()
    function DFJob(name, structure, calculations, local_dir, server, server_dir, header, metadata)
        if !isabspath(local_dir)
            local_dir = abspath(local_dir)
        end
        if isempty(structure.name)
            structure.name = split(name, "_")[1]
        end
        return new(name, structure, calculations, local_dir, server, server_dir, header, metadata)
    end
end

#TODO implement abinit
function DFJob(job_name, structure::AbstractStructure, calculations::Vector{<:DFInput}, common_flags::Pair{Symbol, <:Any}...;
                    pseudoset=:default,
                    pseudospecifier="",
                    job_kwargs...)

    shared_flags = typeof(common_flags) <: Dict ? common_flags : Dict(common_flags...)
    for i in calculations
        i.flags = merge(shared_flags, i.flags)
    end
    out = DFJob(name = job_name, structure = structure, inputs = calculations; job_kwargs...)
    set_pseudos!(out, pseudoset, pseudospecifier, print=false)
    return out
end

DFJob(job_name, ciffile::String, calculations::Vector{<:DFInput}, args...; kwargs...) =
    DFJob(job_name, Structure(ciffile, name = job_name), calculations, args... ; kwargs...)

function DFJob(job::DFJob, flagstoset...; cell_=copy(cell(job)), atoms_=copy(atoms(job)), name=job.name,
                                          server_dir = job.server_dir,
                                          local_dir  = job.local_dir,
                                          pseudoset  = nothing,
                                          pseudospecifier = "")
    newjob = deepcopy(job)

    set_cell!(newjob, cell_)
    if pseudoset === nothing
        pseudoset, specifier = getpseudoset(job.structure.atoms[1])
        specifier = pseudospecifier === nothing ? specifier : pseudospecifier
        set_atoms!(newjob, atoms_, pseudoset = pseudoset, pseudospecifier=specifier)
    else
        set_atoms!(newjob, atoms_, pseudoset = pseudoset, pseudospecifier= pseudospecifier)
    end
    set_serverdir!(newjob, server_dir)
    set_localdir!(newjob, local_dir)
    newjob.name = name

    set_flags!(newjob, flagstoset..., print=false)
    return newjob
end

"""
    DFJob(job_dir::String, job_fuzzy = "job"; kwargs...)

Loads and returns a local DFJob. kwargs will be passed to the constructor.
"""
DFJob(job_dir::String, job_fuzzy="job"; kwargs...) =
    DFJob(;merge(merge((local_dir=job_dir,), read_job_inputs(joinpath(job_dir, searchdir(job_dir, job_fuzzy)[1]))), kwargs)...)

name(job) = job.name
#-------------------BEGINNING GENERAL SECTION-------------#
scriptpath(job::DFJob)       = joinpath(job.local_dir, "job.tt")
starttime(job::DFJob)        = mtime(scriptpath(job))

runslocal(job::DFJob)        = job.server        =="localhost"
structure(job::DFJob)        = job.structure
isQEjob(job::DFJob)          = any(x->package(x) == QE, inputs(job))
iswannierjob(job::DFJob)     = any(x->package(x) == Wannier90, inputs(job)) && any(x->isnscfcalc(x), inputs(job))
getnscfcalc(job::DFJob)      = getfirst(x -> isnscfcalc(x), inputs(job))
cell(job::DFJob)             = cell(structure(job))

input(job::DFJob, n::String) = getfirst(x -> occursin(n, name(x)), inputs(job))
inputs(job::DFJob)           = job.inputs

"""
    inputs(job::DFJob, names::Vector)

Returns an array of the inputs that match the names.
"""
inputs(job::DFJob, names::Vector, fuzzy=true) = fuzzy ? filter(x -> any(occursin.(names, name(x))), inputs(job)) : input.(job, names)
inputs(job::DFJob, n::String, fuzzy=true) = inputs(job, [n], fuzzy)
inputs(job::DFJob, package_::Package) = filter(x->package(x)==package_, inputs(job))
inpath(job::DFJob, n) = inpath(input(job,n))
outpath(job::DFJob, n) = outpath(input(job,n))

"Runs some checks on the set flags for the inputs in the job, and sets metadata (:prefix, :outdir etc) related flags to the correct ones. It also checks whether flags in the various inputs are allowed and set to the correct types."
function sanitizeflags!(job::DFJob)
    set_flags!(job, :prefix => "$(job.name)", print=false)
    if iswannierjob(job)
	    nscfcalc = getnscfcalc(job)
	    if package(nscfcalc) == QE && hasflag(nscfcalc, :nbnd)
	        set_flags!(job, :num_bands => nscfcalc[:nbnd], print=false)
        elseif package(nscfcalc) == Elk
	        setflags!(job, :num_bands => length(nscfcalc[:wann_bands]))
	        nscfcalc[:wann_projections] = projections_string.(unique(filter(x->!isempty(projections(x)), atoms(job))))
	        nscfcalc[:elk2wan_tasks]    = ["602", "604"]
	        nscfcalc[:wann_seedname]    = Symbol(name(job))
	        if job[:wannier_plot] == true
		        push!(nscfcalc[:elk2wan_tasks], "605")
	        end
        end
    end
    for i in filter(x -> package(x) == QE, inputs(job))
        outdir = isempty(job.server_dir) ? joinpath(job, "outputs") : joinpath(job.server_dir, splitdir(job.local_dir)[end], "outputs")
        set_flags!(i, :outdir => "$outdir", print=false)
	    ecutwfc, ecutrho = find_cutoffs(job) # Ideally this should also be at the end stage
	    if exec(i, "pw.x") !== nothing
    	    if !hasflag(i, :ecutwfc)
        	    @info "No energy cutoff was specified in input with name: $(name(i))\nCalculating one from the pseudo files.\nCalculated ecutwfc=$ecutwfc, ecutrho=$ecutrho."
        	    set_flags!(i, :ecutwfc => ecutwfc, :ecuthro => ecutrho, print=false)
    	    end
	    end
    end
    sanitize_pseudos!(job)
    sanitizeflags!.(inputs(job))
end

function sanitize_pseudos!(job::DFJob)
	all_pseudos = pseudo.(atoms(job))
	uni_dirs    = unique(map(x->x.dir, all_pseudos))
	uni_pseudos = unique(all_pseudos)
	pseudo_dir  = length(uni_dirs) == 1 ? uni_dirs[1] : job.local_dir 
	if length(uni_dirs) > 1
		@info "Found pseudos in multiple directories, copying them to job directory"
		for pseudo in uni_pseudos
			cp(path(pseudo), joinpath(job.local_dir, pseudo.name), force=true)
		end
	end
	for p in all_pseudos
		p.dir = pseudo_dir
	end
end

function sanitize_magnetization!(job::DFJob)
    if !any(x->package(x) == QE, inputs(job))
        return
    end
    sanitize_magnetization!(job.structure)
end

function find_cutoffs(job::DFJob)
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
    return maxecutwfc, maxecutrho
end

function sanitize_projections!(job::DFJob)
    if !any(x->!isempty(projections(x)), atoms(job))
        return
    end
    uats = unique(atoms(job))
    projs = unique([name(at) => [p.orb.name for p in projections(at)] for at in uats])
    set_projections!(job, projs...;print=false)
end

"Checks the last created output file for a certain job."
function runninginput(job::DFJob)
    @assert job.server == "localhost" "Intended use for now is locally."
    t = mtime(scriptpath(job))
    for i in reverse(inputs(job))
        p = outpath(i)
        if ispath(p) && mtime(p) > t
            return i
        end
    end
end

"Finds the input corresponding to the name and returns the full output path."
outpath(job::DFJob, n::String) = outpath(input(job,n))

Base.joinpath(job::DFJob, n::AbstractString) = joinpath(job.local_dir, n)

runslocal_assert(job::DFJob) =
    @assert runslocal(job) "This only works if the job runs on `localhost`."

function Base.pop!(job::DFJob, name::String)
    i = findfirst(x -> x.name == name, job.inputs)
    if i === nothing
        error("No input with name $name found.")
    end
    out = job.inputs[i]
    deleteat!(job.inputs, i)
    return out
end

searchdir(job::DFJob, str::AbstractString) = joinpath.((job,), searchdir(job.local_dir, str))

