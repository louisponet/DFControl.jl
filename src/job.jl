#here all the different input structures for the different calculations go
"""
Represents a full DFT job with multiple input files and calculations.
"""
mutable struct DFJob
    id           ::Int
    name         ::String
    structure    ::AbstractStructure
    inputs       ::Vector{DFInput}
    local_dir    ::String
    server       ::String
    server_dir   ::String
    header       ::Vector{String}
    metadata     ::Dict
    function DFJob(name, structure, calculations, local_dir, server, server_dir, header = getdefault_jobheader())
        if local_dir != ""
            local_dir = local_dir
        end

        if server_dir != ""
            server_dir = server_dir
        end
        if !isabspath(local_dir)
            local_dir = abspath(local_dir)
        end
        test = filter(x -> x.name == name,UNDO_JOBS)
        if length(test) == 1
            job = new(test[1].id, name, structure, calculations, local_dir, server, server_dir, header, Dict())
            UNDO_JOBS[test[1].id] = deepcopy(job)
        elseif length(test) == 0
            job = new(length(UNDO_JOBS) + 1, name, structure, calculations, local_dir, server, server_dir, header, Dict())
            push!(UNDO_JOBS, deepcopy(job))
        end
        job
    end
end

#TODO implement abinit
function DFJob(job_name, local_dir, structure::AbstractStructure, calculations::Vector, common_flags...;
                    server=getdefault_server(),
                    server_dir="",
                    package=QE,
                    bin_dir=joinpath("usr","local","bin"),
                    pseudoset=:default,
                    pseudospecifier="",
                    header=getdefault_jobheader())

    @assert package==QE "Only implemented for Quantum Espresso!"

    job_calcs = DFInput[]
    if typeof(common_flags) != Dict
        common_flags = Dict(common_flags)
    end
    bin_dir = bin_dir
    req_flags = Dict(:prefix  => "$job_name",
                     :outdir => "$server_dir",
                     :ecutwfc => 25.)
    merge!(req_flags, common_flags)
    for (calc, (excs, data)) in calculations
        calc_ = typeof(calc) == String ? Symbol(calc) : calc
        if in(calc_, [:vc_relax, :relax, :scf])
            k_points = get(data, :k_points, [1, 1, 1, 0, 0, 0])
            k_option = :automatic
        elseif calc_ == :nscf
            k_points = kgrid(get(data, :k_points, [1, 1, 1])[1:3]..., QE)
            k_option = :crystal
        elseif calc_ == :bands
            k_points = get(data, :k_points, [[0., 0., 0., 1.]])
            num_k = 0.0
            for point in k_points
                num_k += point[4]
            end
            if num_k > 100.
                if !haskey(data, :flags)
                    data[:flags] = Pair{Symbol, Any}[]
                end
                push!(data[:flags], :verbosity => "high")
            end
            k_option = :crystal_b
        end
        flags  = convert(Vector{Pair{Symbol, Any}}, get(data, :flags, Pair{Symbol, Any}[]))
        if excs[2].exec == "pw.x"
            push!(flags, :calculation => "$(string(calc_))")
            datablocks = [InputData(:k_points, k_option, k_points)]
        else
            datablocks =  InputData[]
        end
        input_ = DFInput{package}(string(calc_), local_dir,
                         Dict{Symbol, Any}(),
                         datablocks, excs, true)
        setflags!(input_, req_flags..., print=false)
        setflags!(input_, flags..., print=false)
        push!(job_calcs, input_)
    end
    out = DFJob(job_name, structure, job_calcs, local_dir, server, server_dir, header)
    setatoms!(out, structure.atoms, pseudoset = pseudoset, pseudospecifier= pseudospecifier)
    return DFJob(job_name, structure, job_calcs, local_dir, server, server_dir, header)
end

function DFJob(job_name, local_dir, ciffile::String, calculations::Vector, args...; kwargs...)
    structure = Structure(ciffile, name=job_name)
    return DFJob(job_name, local_dir, structure, calculations, args... ; kwargs...)
end

function DFJob(job::DFJob, flagstoset...; cell_=copy(cell(job)), atoms_=copy(atoms(job)), name=job.name,
                                          server_dir = job.server_dir,
                                          local_dir  = job.local_dir,
                                          pseudoset  = nothing,
                                          pseudospecifier = "")
    newjob = deepcopy(job)

    setcell!(newjob, cell_)
    if pseudoset == nothing
        pseudoset, specifier = getpseudoset(job.structure.atoms[1])
        specifier = pseudospecifier == nothing ? specifier : pseudospecifier
        setatoms!(newjob, atoms_, pseudoset = pseudoset, pseudospecifier=specifier)
    else
        setatoms!(newjob, atoms_, pseudoset = pseudoset, pseudospecifier= pseudospecifier)
    end
    setserverdir!(newjob, server_dir)
    setlocaldir!(newjob, local_dir)
    newjob.name = name

    setflags!(newjob, flagstoset..., print=false)
    return newjob
end

"""
    DFJob(job_dir::String, T=Float64; job_fuzzy = "job", new_job_name=nothing, new_local_dir=nothing, server=getdefault_server(),server_dir="")

Loads and returns a local DFJob. If local_dir is not specified the job directory will be registered as the local one.
"""
function DFJob(job_dir::String, T=Float64;
                  job_fuzzy     = "job",
                  new_job_name  = "",
                  new_local_dir = nothing,
                  server        = getdefault_server(),
                  server_dir    = "")
    name, header, inputs, structure = read_job_inputs(joinpath(job_dir, searchdir(job_dir, job_fuzzy)[1]))
    j_name = isempty(new_job_name) ? name : new_job_name
    structure_name = split(j_name, "_")[1]
    structure.name = structure_name

    if new_local_dir != nothing
        return DFJob(j_name, structure, inputs, new_local_dir, server, server_dir, header)
    else
        return DFJob(j_name, structure, inputs, job_dir, server, server_dir, header)
    end
end

"""
    DFJob(server_dir::String, local_dir::String, server=getdefault_server(); job_fuzzy="*job*", new_job_name="")

Pulls a server job to local directory and then loads it. A fuzzy search for the job file will be performed and the found input files will be pulled.
"""
function DFJob(server_dir::String, local_dir::String, server = getdefault_server();
                         job_fuzzy    = "*job*",
                         new_job_name = "")

    pulljob(server, server_dir, local_dir)
    return DFJob(local_dir, server=server, server_dir=server_dir, new_job_name=new_job_name)
end

#-------------------BEGINNING GENERAL SECTION-------------#
scriptpath(job::DFJob) = joinpath(job.local_dir, "job.tt")
starttime(job::DFJob) = mtime(scriptpath(job))

runslocal(job::DFJob) = job.server=="localhost"
structure(job::DFJob) = job.structure
iswannierjob(job::DFJob) = any(x->package(x) == Wannier90, inputs(job)) && any(x->flag(x, :calculation) == "nscf", inputs(job))
getnscfcalc(job::DFJob) = getfirst(x->flag(x, :calculation) == "nscf", inputs(job))
cell(job::DFJob) = cell(structure(job))

input(job::DFJob, n::String) = getfirst(x -> occursin(n, name(x)), inputs(job))
inputs(job::DFJob) = job.inputs

"""
    inputs(job::DFJob, names::Vector)

Returns an array of the inputs that match the names.
"""
inputs(job::DFJob, names::Vector, fuzzy=true) = fuzzy ? filter(x -> any(occursin.(names, name(x))), inputs(job)) : input.(job, names)
inputs(job::DFJob, n::String, fuzzy=true) = inputs(job, [n], fuzzy)
inputs(job::DFJob, package_::Package) = filter(x->package(x)==package_, inputs(job))

function Base.getindex(job::DFJob, id::String)
    tmp = getfirst(x -> name(x)==id, inputs(job))
    if tmp != nothing
        return tmp
    else
        error("No Input with name $id")
    end
end
function Base.setindex!(job::DFJob, dat, key::Symbol)
    for input in inputs(job)
        input[key] = dat
    end
end

setname!(job::DFJob, oldn, newn) = (input(job, oldn).name = newn)
inpath(job::DFJob, n) = inpath(input(job,n))
outpath(job::DFJob, n) = outpath(input(job,n))

#TODO should we also create a config file for each job with stuff like server etc? and other config things,
#      which if not supplied could contain the default stuff?
"""
    pulljob(server::String, server_dir::String, local_dir::String; job_fuzzy="*job*")

Pulls job from server. If no specific inputs are supplied it pulls all .in and .tt files.
"""
function pulljob(server::String, server_dir::String, local_dir::String; job_fuzzy="*job*")
    server_dir = server_dir
    local_dir  = local_dir
    if !ispath(local_dir)
        mkpath(local_dir)
    end

    pull_server_file(filename) = pullfile(server, server_dir, local_dir, filename)
    pull_server_file(job_fuzzy)
    job_file = searchdir(local_dir, strip(job_fuzzy, '*'))[1]

    if job_file != nothing
        input_files, output_files = read_job_filenames(joinpath(local_dir, job_file))
        for file in input_files
            pull_server_file(file)
        end
    end
end

pulljob(args...; kwargs...) = pulljob(getdefault_server(), args..., kwargs...)


"""
    save(job::DFJob)

Saves a DFJob, it's job file and all it's input files.
"""
function save(job::DFJob, local_dir=job.local_dir)
    local_dir = local_dir != "" ? local_dir : error("Please specify a valid local_dir!")
    if !ispath(local_dir)
        mkpath(local_dir)
        @info "$local_dir did not exist, it was created."
    end
    sanitizeflags!(job)
    job.local_dir = local_dir
    return writejobfiles(job)
end

"Runs some checks on the set flags for the inputs in the job, and sets metadata (:prefix, :outdir etc) related flags to the correct ones. It also checks whether flags in the various inputs are allowed and set to the correct types."
function sanitizeflags!(job::DFJob)
    setflags!(job, :prefix => "$(job.name)", print=false)
    if iswannierjob(job)
        setflags!(job, :num_bands => flag(getnscfcalc(job), :nbnd), print=false)
    end
    sanitizeflags!.(inputs(job))
end

#TODO only uses qsub for now. how to make it more general?
"""
    submit(job::DFJob; server=job.server, server_dir=job.server_dir)

Saves the job locally, and then either runs it locally using `qsub` (when `job.server == "localhost"`) or sends it to the specified `job.server` in `job.server_dir`, and submits it using `qsub` on the server.
"""
function submit(job::DFJob; server=job.server, server_dir=job.server_dir)
    save(job)
    job.server = server
    job.metadata[:slurmid] = qsub(job)
end

"Checks the last created output file for a certain job. Intended use for now is locally."
function runninginput(job::DFJob)
    t = mtime(scriptpath(job))
    for i in reverse(inputs(job))
        p = outpath(i)
        if ispath(p) && mtime(p) > t
            return i
        end
    end
end

"""
    abort(job::DFJob)

Will look for the job id inside it's metadata and try to remove it from the server queue. If the lastrunning input happened to be a QE input, the correct abort file will be written. If it's Wannier90 the job will be brutally removed from the slurm queue. EDIT: It's absolutely impossible to gracefully abort a multi job script with QE... for later
"""
function abort(job::DFJob)
    lastrunning = runninginput(job)
    if lastrunning == nothing
        error("Is this job running?")
    end
    if package(lastrunning) == QE
        writeabortfile(job, lastrunning)
    else
        if !haskey(job.metadata, :slurmid)
            error("No slurm id found for this job.")
        end
        qdel(job::DFJob)
    end
end

"""
    add!(job::DFJob, input::DFInput, index::Int=length(job.inputs)+1; n=name(input))

Adds a calculation to the job, at the specified index.
"""
function add!(job::DFJob, input::DFInput, index::Int=length(job.inputs)+1; n=name(input))
    input.name = n
    insert!(job.inputs, index, input)
    return job
end

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

"""
    flag(job::DFJob, inputs::Vector{<:DFInput}, flag_name::Symbol)

Looks through the input names and returns the value of the specified flag.
"""
function flag(job::DFJob, inputs::Vector{<:DFInput}, fl::Symbol)
    for calc in inputs
        flag_ = flag(calc, fl)
        if flag_ != nothing
            return flag_
        end
    end
    @warn "Flag $fl not found in any input."
end

"""
    flag(job::DFJob, flag_name::Symbol)

Looks through all the calculations and returns the value of the specified flag.
"""
flag(job::DFJob, flag_name::Symbol) =
    flag(job, inputs(job), flag_name)
flag(job::DFJob, name::String, flag_name::Symbol) =
    flag(job, [input(job, name)], flag_name)

"""
    data(job::DFJob, name::String, dataname::Symbol)

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
    setexecflags!(job::DFJob, exec, flags...)

Goes through the calculations of the job and if the name contains any of the `inputnames` it sets the exec flags to the specified ones.
"""
setexecflags!(job::DFJob, exec, flags...) =
    setexecflags!.(job.inputs, (exec, flags)...)
rmexecflags!(job::DFJob, exec, flags...) =
    rmexecflags!.(job.inputs, (exec, flags)...)

"Returns the executables attached to a given input."
execs(job::DFJob, name) = execs(input(job, name))

setexecdir!(job::DFJob, exec, dir) = setexecdir!.(job.inputs, exec, dir)


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
"""
    atoms(job::DFJob)

Returns a list the atoms in the structure.
"""
atoms(job::DFJob) = atoms(job.structure)

"Returns the ith atom with id `atsym`."
atom(job::DFJob, atsym::Symbol, i=1) = filter(x -> x.id == atsym, atoms(job))[i]

#automatically sets the cell parameters for the entire job, implement others
"""
    setcell_parameters!(job::DFJob, cell_::Mat3)

sets the cell parameters of the structure in the job.
"""
function setcell!(job::DFJob, cell_::Mat3)
    job.structure.cell = cell_
    return job
end

"""
    setkpoints!(job::DFJob, n, k_points)

sets the data in the k point `DataBlock` inside the specified inputs.
"""
function setkpoints!(job::DFJob, n, k_points; print=true)
    for calc in inputs(job, n)
        setkpoints!(calc, k_points, print=print)
    end
    return job
end


"sets the pseudopotentials to the specified one in the default pseudoset."
function setpseudos!(job::DFJob, set, specifier="")
    setpseudos!(job.structure, set, specifier)
    dir = getdefault_pseudodir(set)
    dir != nothing && setflags!(job, :pseudo_dir => "$dir", print=false)
    return job
end

"sets the pseudopotentials to the specified one in the default pseudoset."
function setpseudos!(job::DFJob, pseudodir, at_pseudos::Pair{Symbol, String}...)
    setpseudos!(job.structure, at_pseudos...)
    setflags!(job, :pseudo_dir => "$pseudodir", print=false)
    return job
end
"""
    setheaderword!(job::DFJob, word::String, new_word::String)


Replaces the specified word in the header with the new word.
"""
function setheaderword!(job::DFJob, word::String, new_word::String; print=true)
    for (i, line) in enumerate(job.header)
        if occursin(word, line)
            job.header[i] = replace(line, word => new_word)
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


"""
sets the projections of the specified atoms inside the job structure.
"""
setprojections!(job::DFJob, projections...) =
    setprojections!(job.structure, projections...)

"Returns the projections inside the job for the specified `i`th atom in the job with id `atsym`."
projections(job::DFJob, atsym::Symbol, i=1) = projections(atom(job, atsym, i))
"Returns all the projections inside the job."
projections(job::DFJob) = projections.(atoms(job))


"""
    addwancalc!(job::DFJob, nscf::DFInput{QE}, Emin::Real, projections;
                     Emin=-5.0,
                     Epad=5.0,
                     wanflags=SymAnyDict(),
                     pw2wanexec=Exec("pw2wannier90.x", nscf.execs[2].dir, nscf.execs[2].flags),
                     wanexec=Exec("wannier90.x", nscf.execs[2].dir),
                     bands=readbands(nscf))

Adds a wannier calculation to a job. For now only works with QE.
"""
function addwancalc!(job::DFJob, nscf::DFInput{QE}, Emin::Real, projections_...;
                     Epad=5.0,
                     wanflags=SymAnyDict(),
                     pw2wanexec=Exec("pw2wannier90.x", nscf.execs[2].dir),
                     wanexec=Exec("wannier90.x", nscf.execs[2].dir),
                     bands=readbands(nscf),
                     print=true)

    spin = isspincalc(nscf)
    if spin
        pw2wannames = ["pw2wan_up", "pw2wan_dn"]
        wannames = ["wanup", "wandn"]
        print && (@info "Spin polarized calculation found (inferred from nscf input).")
    else
        pw2wannames = ["pw2wan"]
        wannames = ["wan"]
    end

    @assert flag(nscf, :calculation) == "nscf" error("Please provide a valid 'nscf' calculation.")
    if flag(nscf, :nosym) != true
        print && (@info "'nosym' flag was not set in the nscf calculation.\nIf this was not intended please set it and rerun the nscf calculation.\nThis generally gives errors because of omitted kpoints, needed for pw2wannier90.x")
    end

    setprojections!(job, projections_...)
    nbnd = nprojections(job.structure)
    print && (@info "num_bands=$nbnd (inferred from provided projections).")

    wanflags = SymAnyDict(wanflags)
    wanflags[:dis_win_min], wanflags[:dis_froz_min], wanflags[:dis_froz_max], wanflags[:dis_win_max] = wanenergyranges(Emin, nbnd, bands, Epad)

    wanflags[:num_bands] = length(bands)
    wanflags[:num_wann]  = nbnd
    kpoints = data(nscf, :k_points).data
    wanflags[:mp_grid] = kakbkc(kpoints)
    wanflags[:preprocess] = true
    print && (@info "mp_grid=$(join(wanflags[:mp_grid]," ")) (inferred from nscf input).")

    kdata = InputData(:kpoints, :none, [k[1:3] for k in kpoints])

    for (pw2wanfil, wanfil) in zip(pw2wannames, wannames)
        add!(job, DFInput{Wannier90}(wanfil, job.local_dir, copy(wanflags), [kdata], [Exec(), wanexec], true))
    end

    setfls!(job, name, flags...) = setflags!(job, name, flags..., print=false)
    if spin
        setfls!(job, "wanup", :spin => "up")
        setfls!(job, "wandn", :spin => "down")
    end
    return job
end

"""
    addwancalc!(job::DFJob, nscf::DFInput, projwfc::DFInput, threshold::Real, projections...; kwargs...)

Adds a wannier calculation to the job, but instead of passing Emin manually, the output of a projwfc.x run
can be used together with a `threshold` to determine the minimum energy such that the contribution of the
projections to the DOS is above the `threshold`.
"""
function addwancalc!(job::DFJob, nscf::DFInput, projwfc::DFInput, threshold::Real, projections::Pair...; kwargs...)
    @assert hasoutfile(projwfc) @error "Please provide a projwfc Input that has an output file."
    Emin = Emin_from_projwfc(job, outpath(projwfc), threshold, projections...)
    addwancalc!(job, nscf, Emin, projections...; kwargs...)
end

addwancalc!(job::DFJob, nscf_name::String, Emin::Real, projections::Pair...; kwargs...) =
    addwancalc!(job, input(job, nscf_name), Emin, projections...; kwargs...)

addwancalc!(job::DFJob, nscf_name::String, projwfc_name::String, threshold::Real, projections::Pair...; kwargs...) =
    addwancalc!(job, input(job, nscf_name), input(job, projwfc_name), threshold, projections...; kwargs...)


"Automatically calculates and sets the wannier energies. This uses the projections, `Emin` and the bands to infer the other limits.\n`Epad` allows one to specify the padding around the inner and outer energy windows"
function setwanenergies!(job::DFJob, bands, Emin::Real; Epad=5.0, print=true)
    wancalcs = filter(x -> package(x) == Wannier90, job.inputs)
    @assert length(wancalcs) != 0 error("Job ($(job.name)) has no Wannier90 calculations, nothing todo.")
    nbnd = sum([sum(orbsize.(t)) for  t in projections(job)])
    print && (@info "num_bands=$nbnd (inferred from provided projections).")
    winmin, frozmin, frozmax, winmax = wanenergyranges(Emin, nbnd, bands, Epad)
    map(x->setflags!(x, :dis_win_min => winmin, :dis_froz_min => frozmin, :dis_froz_max => frozmax, :dis_win_max => winmax, :num_wann => nbnd, :num_bands=>length(bands); print=false), wancalcs)
    return job
end

function Emin_from_projwfc(job::DFJob, projwfc::String, threshold::Number, projections::Pair...)
    states, bands = read_qe_projwfc(projwfc)
    mask = zeros(length(states))
    for (atsym, projs) in projections
        atids = findall(x -> x.id == atsym, atoms(job))
        stateids = Int[]
        for proj in projs
            orb = orbital(proj)
            push!.((stateids,), findall(x -> x.atom_id ∈ atids && x.l == orb.l, states))
        end
        mask[stateids] .= 1.0
    end
    Emin = 0.0
    for b in bands
        ψ = mean(b.extra[:ψ])
        tot_relevant_occupation = dot(mask, ψ)
        if tot_relevant_occupation > threshold
            Emin = minimum(b.eigvals)
            break
        end
    end
    return Emin
end

"Creates a new `DFInput` from the template with the new flags and new data, then adds it to the inputs of the job at the specified index."
addcalc!(job::DFJob, input::DFInput, index::Int=length(job.inputs)+1) = insert!(job.inputs, index, input)

function addcalc!(job::DFJob, template::DFInput, name::String, newflags...; index=length(job.inputs)+1, run=true, newdata=nothing)
    newcalc = DFInput(template, name, newflags..., data=newdata, run=run)
    addcalc!(job, newcalc, index)
    job
end

"""
    addcalc!(job::DFJob, template::DFInput, kpoints::Vector{NTuple{4}}, newflags...; name="bands", run=true, template="scf")

Searches for the given template and creates a bands calculation from it.
"""
function addcalc!(job::DFJob, template::DFInput, kpoints::Vector{<:NTuple{4}}, args...; name="bands", kwargs...)
    addcalc!(job, template, name, :calculation => "bands",args...; kwargs...)
    setkpoints!(job, name, kpoints, print=false)
    job
end

"""
    addcalc!(job::DFJob, template::DFInput, kpoints::NTuple{3}, newflags...; name="nscf", run=true, template="scf")

Searches for the given template and creates a bands calculation from it.
"""
function addcalc!(job::DFJob, template::DFInput, kpoints::NTuple{3}, args...; name="nscf", kwargs...)
    addcalc!(job, template, name, :calculation => "nscf",args...; kwargs...)
    setkpoints!(job, name, kpoints, print=false)
    job
end
"""
    addcalc!(job::DFJob, template::DFInput, kpoints::NTuple{6}, newflags...; name="scf", run=true, template="nscf")

Searches for the given template and creates a bands calculation from it.
"""
function addcalc!(job::DFJob, template::DFInput, kpoints::NTuple{6}, args...; name="scf", kwargs...)
    addcalc!(job, template, name, :calculation => "scf", args...; kwargs...)
    setkpoints!(job, name, kpoints, print=false)
    job
end

addcalc!(job::DFJob, template::String, args...; kwargs...) = addcalc!(job, input(job, template), args...; kwargs...)
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

"Finds the input corresponding to the name and returns the full output path."
outpath(job::DFJob, n::String) = outpath(input(job,n))

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


function isrunning(job::DFJob)
    @assert haskey(job.metadata, :slurmid) error("No slurmid found for job $(job.name)")
    cmd = `qstat -f $(job.metadata[:slurmid])`
    if runslocal(job)
        str = read(cmd, String)
    else
        str = sshreadstring(job.server, cmd)
    end
    isempty(str) && return false
    splstr = split(str)
    for (i,s) in enumerate(splstr)
        if s=="job_state"
            return any(splstr[i+2] .== ["Q","R"])
        end
    end
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

"Reads throught the pseudo files and tries to figure out the correct cutoffs"
function setcutoffs!(job::DFJob)
    @assert job.server == "localhost" "Cutoffs can only be automatically set if the pseudo files live on the local machine."
    pseudofiles = filter(!isempty, [pseudo(at) for at in atoms(job)])
    pseudodirs  = String[]
    for i in inputs(job)
        if package(i) == QE
            dr = pseudodir(i)
            if dr != nothing && ispath(dr) #absolute paths only allowed in QE
                push!(pseudodirs, dr)
            end
        end
    end
    @assert !isempty(pseudofiles) "No atoms with pseudo files found."
    @assert !isempty(pseudodirs) "No valid pseudo directories found in the inputs."
    maxecutwfc = 0.0
    maxecutrho = 0.0
    for d in pseudodirs
        for f in pseudofiles
            pth = joinpath(d, f)
            if ispath(pth)
                println(pth)
                ecutwfc, ecutrho = read_cutoffs_from_pseudofile(pth)
                show(ecutwfc)
                show(ecutrho)
                if ecutwfc != nothing && ecutrho != nothing
                    maxecutwfc = ecutwfc > maxecutwfc ? ecutwfc : maxecutwfc
                    maxecutrho = ecutrho > maxecutrho ? ecutrho : maxecutrho
                end
            end
        end
    end
    setcutoffs!.(inputs(job), maxecutwfc, maxecutrho)
end
