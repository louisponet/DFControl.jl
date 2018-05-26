#here all the different input structures for the different calculations go
"""
Represents a full DFT job with multiple input files and calculations.
"""
mutable struct DFJob
    id           ::Int
    name         ::String
    structure    ::AbstractStructure
    calculations ::Vector{DFInput}
    local_dir    ::String
    server       ::String
    server_dir   ::String
    header       ::Vector{String}
    function DFJob(name, structure, calculations, local_dir, server, server_dir, header = getdefault_jobheader())
        if local_dir != ""
            local_dir = form_directory(local_dir)
        end

        if server_dir != ""
            server_dir = form_directory(server_dir)
        end

        test = filter(x -> x.name == name,UNDO_JOBS)
        if length(test) == 1
            job = new(test[1].id, name, structure, calculations, local_dir, server, server_dir, header)
            UNDO_JOBS[test[1].id] = deepcopy(job)
        elseif length(test) == 0
            job = new(length(UNDO_JOBS) + 1, name, structure, calculations, local_dir, server, server_dir, header)
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
                    bin_dir="/usr/local/bin/",
                    pseudo_set=:default,
                    pseudo_specifier="",
                    header=getdefault_jobheader())

    @assert package==QE "Only implemented for Quantum Espresso!"
    local_dir = form_directory(local_dir)
    job_calcs = DFInput[]
    if typeof(common_flags) != Dict
        common_flags = Dict(common_flags)
    end
    bin_dir = form_directory(bin_dir)
    req_flags = Dict(:prefix  => "'$job_name'",
                     :outdir => "'$server_dir'",
                     :ecutwfc => 25.)
    merge!(req_flags, common_flags)
    for (calc, (execs, data)) in calculations
        calc_ = typeof(calc) == String ? Symbol(calc) : calc
        if in(calc_, [Symbol("vc-relax"), :relax, :scf])
            k_points = get(data, :k_points, [1, 1, 1, 0, 0, 0])
            k_option = :automatic
        elseif calc_ == :nscf
            k_points = kgrid(get(data, :k_points, [1, 1, 1])..., :nscf)
            k_option = :crystal
        elseif calc_ == :bands
            k_points = get(data, :k_points, [[0., 0., 0., 1.]])
            num_k = 0.0
            for point in k_points
                num_k += point[4]
            end
            if num_k > 100.
                push!(data[:flags], :verbosity => "'high'")
            end
            k_option = :crystal_b
        end
        flags  = get(data, :flags, Dict{Symbol, Any}())
        if execs[2].exec == "pw.x"
            push!(flags, :calculation => "'$(string(calc_))'")
            datablocks = [InputData(:k_points, k_option, k_points)]
        else
            datablocks =  InputData[]
        end
        input_ = DFInput{package}(string(calc_) * ".in",
                         Dict{Symbol, Any}(),
                         datablocks, execs, true)
        setflags!(input_, req_flags..., print=false)
        setflags!(input_, flags..., print=false)
        push!(job_calcs, input_)
    end
    out = DFJob(job_name, structure, job_calcs, local_dir, server, server_dir, header)
    setatoms!(out, structure.atoms, pseudo_set = pseudo_set, pseudo_specifier= pseudo_specifier)
    return DFJob(job_name, structure, job_calcs, local_dir, server, server_dir, header)
end

function DFJob(job_name, local_dir, ciffile::String, calculations::Vector, args...; kwargs...)
    structure = Structure(ciffile, name=job_name)
    return DFJob(job_name, local_dir, structure, calculations, args... ; kwargs...)
end

function DFJob(job::DFJob, flagstoset...;
               cell             = nothing,
               atoms            = nothing,
               pseudo_set       = nothing,
               pseudo_specifier = "",
               server_dir       = nothing,
               local_dir        = nothing,
               name             = nothing)

    newjob = deepcopy(job)

    if cell != nothing
        setcell!(newjob, cell)
    end
    if atoms != nothing
        if pseudo_set == nothing
            pseudo_set, specifier = getpseudoset(job.structure.atoms[1])
            specifier = pseudo_specifier == nothing ? specifier : pseudo_specifier
            setatoms!(newjob, atoms, pseudo_set = pseudo_set, pseudo_specifier=specifier)
        else
            setatoms!(newjob, atoms, pseudo_set = pseudo_set, pseudo_specifier= pseudo_specifier)
        end
    end
    if pseudo_set != nothing && atoms == nothing
        setpseudos!(newjob, pseudo_set, pseudo_specifier)
    end
    if server_dir != nothing
        setserverdir!(newjob, server_dir)
    end
    if local_dir != nothing
        setlocaldir!(newjob, local_dir)
    end
    if name != nothing
        newjob.name = name
    end

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

    job_dir = form_directory(job_dir)

    name, header, inputs, structure = read_job_inputs(job_dir * searchdir(job_dir, job_fuzzy)[1])
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
#all inputs return arrays, input returns the first element if multiple are found
"""
    inputs(job::DFJob, filenames::Array)

Returns an array of the inputs that match one of the filenames.
"""
function inputs(job::DFJob, filenames::Vector)
    out = DFInput[]
    for name in filenames
        push!(out, filter(x -> contains(x.filename, name), job.calculations)...)
    end
    return out
end

"""
    inputs(job::DFJob, filename::String)

Returns an array of the input that matches the filename.
"""
function inputs(job::DFJob, filename::String)
    return filter(x -> contains(x.filename, filename), job.calculations)
end

"""
    input(job::DFJob, filename::String)

Returns the input that matches the filename.
"""
input(job::DFJob, filename::String) = getfirst(x -> contains(x.filename, filename), job.calculations)

"""
    input(job::DFJob, filenames::Array)

Returns an array of the inputs that match one of the filenames.
"""
function input(job::DFJob, filenames::Vector{String})
    return inputs(job, filenames)
end

#TODO should we also create a config file for each job with stuff like server etc? and other config things,
#      which if not supplied could contain the default stuff?
"""
    pulljob(server::String, server_dir::String, local_dir::String; job_fuzzy="*job*")

Pulls job from server. If no specific inputs are supplied it pulls all .in and .tt files.
"""
function pulljob(server::String, server_dir::String, local_dir::String; job_fuzzy="*job*")
    server_dir = form_directory(server_dir)
    local_dir  = form_directory(local_dir)
    if !ispath(local_dir)
        mkpath(local_dir)
    end

    pull_server_file(filename) = pull_file(server, server_dir, local_dir, filename)
    pull_server_file(job_fuzzy)
    job_file = searchdir(local_dir, strip(job_fuzzy, '*'))[1]

    if job_file != nothing
        input_files, output_files = read_job_filenames(local_dir * job_file)
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
    local_dir = local_dir != "" ? form_directory(local_dir) : error("Please specify a valid local_dir!")
    if !ispath(local_dir)
        mkpath(local_dir)
    end
    job.local_dir = local_dir
    return write_job_files(job)
end

#Incomplete everything is hardcoded for now still!!!! make it configurable
"""
    push(job::DFJob)

Pushes a DFJob from it's local directory to its server side directory.
"""
function push(job::DFJob, newfiles)
    if !isfile(job.local_dir * "job.tt")
        save(job)
    end
    try
        run(`ssh -t $(job.server) touch $(job.server_dir * "tmp.in")`)
        run(`ssh -t $(job.server) rm $(job.server_dir * "tmp.in")`)
    catch
        run(`ssh -t $(job.server) mkdir $(job.server_dir)`)
        info("$(job.server_dir) did not exist on $(job.server), it was created")
    end
    for file in newfiles
        run(`scp $(job.local_dir * file) $(job.server * ":" * job.server_dir)`)
    end
    run(`scp $(job.local_dir * "job.tt") $(job.server * ":" * job.server_dir)`)
end

#TODO only uses qsub for now. how to make it more general?
"""
    submit(job::DFJob; server=nothing, server_dir=nothing)

Submit a DFJob. First saves it locally, pushes it to the server then runs the job file on the server.
"""
function submit(job::DFJob; server=job.server, server_dir=job.server_dir)
    new_files = save(job)
    if server != ""
        push(job, new_files)
        run(`ssh -t $(job.server) cd $(job.server_dir) '&&' qsub job.tt`)
    else
        try
            run(`cd $(job.local_dir) '&&' qsub job.tt`)
        catch
            error("Tried submitting on the local machine but got an error executing `qsub`.")
        end
    end
end

"""
    add!(job::DFJob, input::DFInput, index::Int=length(job.calculations)+1; filename=input.filename)

Adds a calculation to the job, at the specified index.
"""
function add!(job::DFJob, input::DFInput, index::Int=length(job.calculations)+1; filename = input.filename)

    UNDO_JOBS[job.id] = deepcopy(job)

    input.filename = filename
    insert!(job.calculations, index, input)
    return job
end

"""
    setflags!(job::DFJob, calculations::Vector{String}, flags...; print=true)

Sets the flags in the calculations to the flags specified.
This only happens if the specified flags are valid for the calculations.
If necessary the correct control block will be added to the calculation (e.g. for QEInputs).

The values that are supplied will be checked whether they are valid.
"""
function setflags!(job::DFJob, calculations::Vector{String}, flags...; print=true)
    UNDO_JOBS[job.id] = deepcopy(job)
    found_keys = Symbol[]

    for calc in input.(job, calculations)
        t_found_keys, = setflags!(calc, flags..., print=print)
        for key in t_found_keys
            if !(key in found_keys) push!(found_keys, key) end
        end
    end

    n_found_keys = Symbol[]
    for (k, v) in flags
        if !(k in found_keys) push!(n_found_keys, k) end
    end
    if print
        if 1 < length(n_found_keys)
            dfprintln("flags '$(join(":" .* String.(n_found_keys),", "))' were not found in the allowed input variables of the specified calculations!")
        elseif length(n_found_keys) == 1
            dfprintln("flag '$(":" * String(n_found_keys[1]))' was not found in the allowed input variables of the specified calculations!")
        end
    end
    return job
end
setflags!(job::DFJob, flags...;kwargs...)                   = setflags!(job, [calc.filename for calc in job.calculations], flags...;kwargs...)
setflags!(job::DFJob, filename::String, flags...;kwargs...) = setflags!(job, [filename], flags...;kwargs...)

"""
    flag(job::DFJob, calc_filenames, flag_name::Symbol)

Looks through the calculation filenames and returns the value of the specified flag.
"""
function flag(job::DFJob, calc_filenames, fl::Symbol)
    for calc in inputs(job, calc_filenames)
        flag_ = flag(calc, fl)
        if flag_ != nothing
            return flag_
        end
    end
    error("Flag $flag not found in any input files.")
end

"""
    flag(job::DFJob, flag_name::Symbol)

Looks through all the calculations and returns the value of the specified flag.
"""
function flag(job::DFJob, flag_name::Symbol)
    for calc in job.calculations
        flag_ = flag(calc, flag_name)
        if flag_ != nothing
            return flag_
        end
    end
    error("Flag $flag_name not found in any input files.")
end

#TODO set so calculations also have a name.
#TODO set after implementing k_point set so you don't need to specify all this crap
"""
    data(job::DFJob, calc_filenames, name::Symbol)

Looks through the calculation filenames and returns the data with the specified symbol.
"""
function data(job::DFJob, calc_filenames, name::Symbol)
    for calc in inputs(job, calc_filenames)
        return data(calc, name)
    end
end

"""
    setdata!(job::DFJob, calc_filenames, data_block_name::Symbol, new_block_data; option=nothing)

Looks through the calculation filenames and sets the data of the datablock with `data_block_name` to `new_block_data`.
if option is specified it will set the block option to it.
"""
function setdata!(job::DFJob, calc_filenames, data_block_name::Symbol, new_block_data; option=nothing)
    UNDO_JOBS[job.id] = deepcopy(job)

    for calc in inputs(job, calc_filenames)
        setdata!(calc, data_block_name, new_block_data, option=option)
    end
end

"""
    rmflags!(job::DFJob, calc_filenames, flags...)

Looks through the calculation filenames and removes the specified flags.
"""
function rmflags!(job::DFJob, calc_filenames::Vector{<:AbstractString}, flags...)
    UNDO_JOBS[job.id] = deepcopy(job)

    for calc in inputs(job, calc_filenames)
        rmflags!(calc, flags...)
    end
    return job
end
rmflags!(job::DFJob, filename::String, flags...) = rmflags!(job, [filename], flags...)
rmflags!(job::DFJob, flags...) = rmflags!(job, [calc.filename for calc in job.calculations], flags...)

"""
    setflow!(job::DFJob, should_runs...)

Sets whether or not calculations should be run. Calculations are specified using their indices.
"""
function setflow!(job::DFJob, should_runs...)
    UNDO_JOBS[job.id] = deepcopy(job)

    for (filename, run) in should_runs
        for input in inputs(job, filename)
            input.run = run
        end
    end
    return job
end

setflow!(job::DFJob, should_runs::Vector{Bool}) = setflow!(job, [calc.filename => run for (calc, run) in zip(job.calculations, should_runs)]...)

"""
    setflow!(job::DFJob, filenames::Array{String,1}, should_run)

Goes throug the calculation filenames and sets whether it should run or not.
"""
setflow!(job::DFJob, filenames::Vector{String}, should_run) = setflow!.(inputs(job, filenames), should_run)

"""
    setexecflags!(job::DFJob, inputnames, flags...)

Goes through the calculations of the job and if the name contains any of the `inputnames` it sets the exec flags to the specified ones.
"""
setexecflags!(job::DFJob, inputnames, flags...) = setexecflags!.(inputs(job, inputnames), flags...)

"""
    setdata!(job::DFJob, filenames, block::Block)

Adds a block to the specified filenames.
"""
function setdata!(job::DFJob, filenames, data::InputData)
    UNDO_JOBS[job.id] = deepcopy(job)

    for input in inputs(job, filenames)
        setdata!(input, data)
    end
    return job
end

"""
    setdata!(job::DFJob, filenames, name, data, option=:none)

Adds a block to the specified filenames.
"""
function setdata!(job::DFJob, filenames, name, data, option=:none)
    UNDO_JOBS[job.id] = deepcopy(job)

    for input in inputs(job, filenames)
        setdata!(input, name, data, option)
    end
    return job
end

"""
    setfilename!(job::DFJob, old_filename::String, new_filename::String)

sets the filename from the old to the new one.
"""
function setfilename!(job::DFJob, old_filename::String, new_filename::String)
    input          = input(job, old_filename)
    old_filename   = input.filename
    input.filename = new_filename
    dfprintln("Input filename in job $(job.name) has been setd from '$old_filename' to '$new_filename'.")
end

#---------------------------------END GENERAL SECTION ------------------#

"""
    setatoms!(job::DFJob, atoms::Dict{Symbol,<:Array{<:Point3,1}}, pseudo_setname=:default, pseudo_specifier=nothing, option=:angstrom)

Sets the data data with atomic positions to the new one. This is done for all calculations in the job that have that data.
If default pseudopotentials are defined, a set can be specified, together with a fuzzy that distinguishes between the possible multiple pseudo strings in the pseudo set.
These pseudospotentials are then set in all the calculations that need it.
All flags which specify the number of atoms inside the calculation also gets set to the correct value.
"""
setatoms!(job::DFJob, atoms::Dict{Symbol,<:Vector{<:Point3}}; kwargs...) = setatoms(job, convert_2atoms(atoms); kwargs...)

function setatoms!(job::DFJob, atoms::Vector{<:AbstractAtom}; pseudo_set = :default, pseudo_specifier="")
    UNDO_JOBS[job.id] = deepcopy(job)

    job.structure.atoms = atoms
    setpseudos!(job, pseudo_set, pseudo_specifier)
    return job
end

"""
    atoms(job::DFJob)

Returns a list the atoms in the structure.
"""
atoms(job::DFJob) = job.structure.atoms

"Returns the ith atom with id `atsym`."
atom(job::DFJob, atsym::Symbol, i=1) = filter(x -> x.id == atsym, atoms(job))[i]
cell(job::DFJob)  = job.structure.cell
#automatically sets the cell parameters for the entire job, implement others
"""
    setcell_parameters!(job::DFJob, cell_param::Matrix)

sets the cell parameters of the structure in the job.
"""
function setcell!(job::DFJob, cell_param::Matrix)
    UNDO_JOBS[job.id] = deepcopy(job)

    job.structure.cell = cell_param
    return job
end

"""
    setkpoints!(job::DFJob, calc_filename, k_points)

sets the data in the k point `DataBlock` inside the specified calculation.
"""
function setkpoints!(job::DFJob, calc_filenames, k_points)
    UNDO_JOBS[job.id] = deepcopy(job)
    for calc in inputs(job, calc_filenames)
        setkpoints!(calc, k_points)
    end
    return job
end

"""
    setdataoption!(job::DFJob, filenames::Array{String,1}, name::Symbol, option::Symbol)

sets the option of specified data block in the specified calculations.
"""
function setdataoption!(job::DFJob, filenames::Vector{String}, name::Symbol, option::Symbol)
    UNDO_JOBS[job.id] = deepcopy(job)

    for calc in inputs(job, filenames)
        setdataoption!(calc, name, option)
    end
    return job
end
setdataoption!(job::DFJob, filename::String, name::Symbol, option::Symbol) = setdataoption!(job, [filename], name, option)

"""
    setdataoption!(job::DFJob, name::Symbol, option::Symbol)

sets the option of specified data block in all calculations that have the block.
"""
setdataoption!(job::DFJob, name::Symbol, option::Symbol) = setdataoption!(job, [i.filename for i in job.calculations], name, option)

"sets the pseudopotentials to the specified one in the default pseudo_set."
function setpseudos!(job::DFJob, pseudo_set, pseudo_specifier="")
    for (i, atom) in enumerate(job.structure.atoms)
        pseudo = getdefault_pseudo(id(atom), pseudo_set, pseudo_specifier=pseudo_specifier)
        if pseudo == nothing
            warn("Pseudo for $(id(atom)) at index $i not found in set $pseudo_set.\n Setting pseudo_dir to './', this should be changed before trying to run the job!")
            setflags!(job, :pseudo_dir => "'./'", print=false)
        else
            setpseudo!(atom, pseudo)
            setflags!(job, :pseudo_dir => "'$(getdefault_pseudodir(pseudo_set))'", print=false)
        end
    end
    return job
end

"""
    replace_header_word!(job::DFJob, word::String, new_word::String)


Replaces the specified word in the header with the new word.
"""
function setheaderword!(job::DFJob, word::String, new_word::String)
    UNDO_JOBS[job.id] = deepcopy(job)

    for (i, line) in enumerate(job.header)
        if contains(line, word)
            job.header[i] = replace(line, word, new_word)
            s = """Old line:
            $line
            New line:
            $(job.header[i])
            """
            dfprintln(s)
        end
    end
    return job
end

"""
    errors(job::DFJob)

Prints the possible error messages in outputs of the `DFJob`.
"""
function errors(job::DFJob)
    out = outputs(job)
    errors  = Dict{String, Vector{String}}()
    for o in out
        errors[o] = read_errors(o)
    end

    for (filename, errs) in errors
        dfprintln("Error in output '$filename':")
        for err in errs
            dfprintln("$err")
        end
    end

    if isempty(errors)
        dfprintln("No errors found for job '$(job.name)'.")
    end
end

"Returns the full path of the input"
path(job::DFJob, inp::DFInput) =
    joinpath(job.local_dir, inp.filename)
path(job::DFJob, inp::String) = path(job::DFJob, input(job, inp))

"Provides the path of the output to the given input, provided that it is pulled already"
function outpath(job::DFJob, input::DFInput{QE})
    opath = splitext(path(job, input))[1] * ".out"
    ispath(opath) ? opath : error("No output file found at $opath,\n please pull the outputs of the job first.")
end
outpath(job::DFJob, inp::String) = outpath(job, input(job, inp))

"""
sets the projections of the specified atoms inside the job structure.
"""
setprojections!(job::DFJob, projections...) = setprojections!(job.structure, projections...)

"Returns the projections inside the job for the specified `i`th atom in the job with id `atsym`."
projections(job::DFJob, atsym::Symbol, i=1) = projections(atom(job, atsym, i))
"Returns all the projections inside the job."
projections(job::DFJob) = projections.(atoms(job))

"""
    addwancalc!(job::DFJob, nscf::DFInput{QE}, projections;
                     Emin=-5.0,
                     Epad=5.0,
                     wanflags=SymAnyDict(),
                     pw2wanexec=Exec("pw2wannier90.x", nscf.execs[2].dir, nscf.execs[2].flags),
                     wanexec=Exec("wannier90.x", nscf.execs[2].dir),
                     bands=read_qe_bands_file(outpath(job, nscf)))

Adds a wannier calculation to a job. For now only works with QE.
"""
function addwancalc!(job::DFJob, nscf::DFInput{QE}, projections;
                     Emin=-5.0,
                     Epad=5.0,
                     wanflags=SymAnyDict(),
                     pw2wanexec=Exec("pw2wannier90.x", nscf.execs[2].dir),
                     wanexec=Exec("wannier90.x", nscf.execs[2].dir),
                     bands=read_qe_bands_file(outpath(job, nscf)))


    spin = spincalc(nscf)
    if spin
        pw2wanfiles = ["pw2wan_up.in", "pw2wan_dn.in"]
        wanfiles = ["wan_up.win", "wan_dn.win"]
        info("Spin polarized calculation found (inferred from nscf input).")
    else
        pw2wanfiles = ["pw2wan.in"]
        wanfiles = ["wan.win"]
    end

    @assert flag(nscf, :calculation) == "'nscf'" error("Please provide a valid 'nscf' calculation.")
    if flag(nscf, :nosym) != true
        info("'nosym' flag was not set in the nscf calculation.\nIf this was not intended please set it and rerun the nscf calculation.\nThis generally gives errors because of omitted kpoints, needed for pw2wannier90.x")
    end

    ats = atoms(job)
    projections = Dict(projections)
    ids   = [id(at) for at in ats]
    @assert all(haskey.(projections, unique(ids))) error("Please supply a set of projections for each unique atom: $(join(unids," ")).")


    nbnd = sum(sum.([[orbsize(orb) for orb in projections[id]] for id in ids]))
    info("num_bands=$nbnd (inferred from provided projections).")


    wanflags = SymAnyDict(wanflags)
    wanflags[:dis_win_min], wanflags[:dis_froz_min], wanflags[:dis_froz_max], wanflags[:dis_win_max] = wanenergyranges(Emin, nbnd, bands, Epad)

    wanflags[:num_bands] = length(bands)
    wanflags[:num_wann]  = nbnd

    wanflags[:mp_grid] = kakbkc(data(nscf, :k_points).data)
    info("mp_grid=$(join(wanflags[:mp_grid]," ")) (inferred from nscf input).")

    pw2wanflags = SymAnyDict(:prefix => flag(nscf, :prefix), :outdir => flag(nscf, :outdir) == nothing ? "'./'" : flag(nscf, :outdir))
    if haskey(wanflags, :write_hr)
        pw2wanflags[:write_amn] = true
        pw2wanflags[:write_mmn] = true
    end
    if haskey(wanflags, :wannier_plot)
        pw2wanflags[:write_unk] = true
    end

    kdata = InputData(:kpoints, :none, kgrid(wanflags[:mp_grid]..., :wan))

    for (pw2wanfil, wanfil) in zip(pw2wanfiles, wanfiles)
        add!(job, DFInput{Wannier90}(wanfil, copy(wanflags), [kdata], [Exec(), wanexec], true))
        add!(job, DFInput{QE}(pw2wanfil, copy(pw2wanflags), InputData[], [nscf.execs[1], pw2wanexec], true))
    end

    setfls!(job, name, flags...) = setflags!(job, name, flags..., print=false)
    if spin
        setfls!(job, "pw2wan_up", :spin_component => "'up'")
        setfls!(job, "pw2wan_dn", :spin_component => "'down'")
        setfls!(job, "wan_up.win", :spin => "'up'")
        setfls!(job, "wan_dn.win", :spin => "'down'")
    end
    addprojections!(projections, ats)

    return job
end


"Automatically calculates and sets the wannier energies. This uses the projections, `Emin` and the bands to infer the other limits.\n`Epad` allows one to specify the padding around the inner and outer energy windows"
function setwanenergies!(job::DFJob, Emin, bands; Epad=5.0)
    wancalcs = filter(x -> eltype(x) == Wannier90, job.calculations)
    @assert length(wancalcs) != 0 error("Job ($(job.name)) has no Wannier90 calculations, nothing todo.")
    nbnd = sum([sum(orbsize.(t)) for  t in projections(job)])
    info("num_bands=$nbnd (inferred from provided projections).")
    winmin, frozmin, frozmax, winmax = wanenergyranges(Emin, nbnd, bands, Epad)
    setflags!.(wancalcs, :dis_win_min => winmin, :dis_froz_min => frozmin, :dis_froz_max => frozmax, :dis_win_max => winmax, print=false)
end

"""
    undo!(job::DFJob)

Undos the last set to the calculations of the job.
"""
function undo!(job::DFJob)
    job.calculations[:] = UNDO_JOBS[job.id].calculations[:]
    dfprintln("Restored the calculations of job '$(job.name)' to their previous state.")
    return job
end

"""
    undo(job::DFJob)

Returns the previous state of the job.
"""
function undo(job::DFJob)
    return deepcopy(UNDO_JOBS[job.id])
end

"""
    addbandscalculation!(job::DFJob, k_path::Vector{Vector{T}}; filename="bands.in", run=true) where T<:AbstractFloat

Checks if there is an scf calculation in the job and takes it's inputs to generate a bands calculation along the given k-path.
"""
function addbandscalculation!(job::DFJob, k_path::Vector{Vector{T}}; filename="bands.in", run=true) where T<:AbstractFloat
    calc = getfirst(x -> eltype(x) == QE && flag(x, :calculation) == "'scf'", job.calculations)
    bands_calc = DFInput{QE}(calc, filename, run=run, data=[:k_points => (:crystal_b, k_path)])
    push!(job.calculations, bands_calc)
    return job
end


"""
Sets the server dir of the job.
"""
function setserverdir!(job, dir)
    dir = form_directory(dir)
    job.server_dir = dir
    return job
end

"""
Sets the local dir of the job.
"""
function setlocaldir!(job, dir)
    dir = form_directory(dir)
    job.local_dir = dir
    return job
end


"""
    print_info(job::DFJob, filenames::Vector{String})

Prints general info of the job.
"""
function print_info(job::DFJob, filenames::Vector{String})
    s = """--------------------
    DFJob:      $(job.name)
    Local_dir:  $(job.local_dir)
    Server:     $(job.server)
    Server_dir: $(job.server_dir)
    $(length(job.calculations)) calculations
    --------------------
    """
    dfprintln(s)
end
print_info(job::DFJob)                   = print_info(job, [calc.filename for calc in job.calculations])
print_info(job::DFJob, filename::String) = print_info(job, [filename])
