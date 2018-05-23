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
                    bin_dir="~/bin/",
                    runcommand=Exec("mpirun", bin_dir, Dict(:np => 24)),
                    pseudo_set=:default,
                    pseudo_specifier="",
                    header=getdefault_jobheader())

    @assert package==:qe "Only implemented for Quantum Espresso!"
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
    for (calc, (exec, data)) in calculations
        exec = typeof(exec) == String ? exec : string(exec)
        calc_ = typeof(calc) == String ? Symbol(calc) : calc
        if in(calc_, [Symbol("vc-relax"), :relax, :scf])
            k_points = get(data, :k_points, [1, 1, 1, 0, 0, 0])
            k_option = :automatic
        elseif calc_ == :nscf
            k_points = get(data, :k_points, (1, 1, 1))
            k_grid   = kgrid(k_points..., :nscf)
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
        if exec == "pw.x"
            push!(flags, :calculation => "'$(string(calc_))'")
            datablocks = [InputInfo(:k_points, k_option, k_points)]
        else
            datablocks =  InputInfo[]
        end
        input_ = DFInput{package}(string(calc_) * ".in",
                         Dict{Symbol, Any}(),
                         datablocks,
                         runcommand, Exec(string(exec), bin_dir),
                         true)
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

    setflags!(newjob, flagstoset...)
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

    name, header, inputs, structure = read_job_inputs(job_dir * search_dir(job_dir, job_fuzzy)[1])
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
function input(job::DFJob, filename::String)
    return getfirst(x -> contains(x.filename, filename), job.calculations)
end

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
    job_file = search_dir(local_dir, strip(job_fuzzy, '*'))[1]

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
    add!(job::DFJob, input::DFInput, index::Int=length(job.calculations)+1; runcommand=input.runcommand, filename=input.filename)

Adds a calculation to the job, at the specified index.
"""
function add!(job::DFJob, input::DFInput, index::Int=length(job.calculations)+1;
                          runcommand = input.runcommand,
                          filename    = input.filename)

    UNDO_JOBS[job.id] = deepcopy(job)

    input.filename = filename
    input.runcommand = runcommand
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
    for calc in inputs(job, calculations)
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
setflags!(job::DFJob, flags...)                   = setflags!(job, [calc.filename for calc in job.calculations], flags...)
setflags!(job::DFJob, filename::String, flags...) = setflags!(job, [filename], flags...)

"""
    flag(job::DFJob, calc_filenames, flag_name::Symbol)

Looks through the calculation filenames and returns the value of the specified flag.
"""
function flag(job::DFJob, calc_filenames, flag_name::Symbol)
    for calc in inputs(job, calc_filenames)
        flag_ = flag(calc, flag_name)
        if flag_ != nothing
            return flag_
        end
    end
    error("Flag $flag_name not found in any input files.")
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
    inputinfo(job::DFJob, calc_filenames, name::Symbol)

Looks through the calculation filenames and returns the block with the specified symbol.
"""
function inputinfo(job::DFJob, calc_filenames, name::Symbol)
    for calc in inputs(job, calc_filenames)
        return inputinfo(calc, name)
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
    setruncommand!(job::DFJob, inputnames, runcommand::Exec)

Goes through the job calculations and if it contains one of the inputnames it sets the run command of the calculation.
"""
function setruncommand!(job::DFJob, inputnames, runcommand::Exec)
    UNDO_JOBS[job.id] = deepcopy(job)

    for calc in inputs(job, inputnames)
        calc.runcommand = runcommand
        dfprintln("Run command of file '$(calc.filename)' is now: '$(calc.runcommand)'")
    end

    return job
end

"""
    runcommand(job::DFJob, inputname)

Returns the `runcommand`.
"""
function runcommand(job::DFJob, inputname)
    for calc in inputs(job, inputname)
        return calc.runcommand
    end
end

"""
    setrunflags!(job::DFJob, inputnames, flags...)

Goes through the calculations of the job and if the name contains any of the `inputnames` it sets the runcommand flags to the specified ones.
"""
function setrunflags!(job::DFJob, inputnames, flags...)
    calcs = inputs(job, inputnames)
    for calc in calcs
        for (f,v) in flags
            calc.runcommand.flags[f] = v
        end
        dfprintln("run flags of calculation $(calc.filename) are now $(calc.runcommand.flags).")
    end
end

"""
    setexecflags!(job::DFJob, inputnames, flags...)

Goes through the calculations of the job and if the name contains any of the `inputnames` it sets the exec flags to the specified ones.
"""
setexecflags!(job::DFJob, inputnames, flags...) = setexecflags!.(inputs(job, inputnames), flags...)

"""
    setinputinfo!(job::DFJob, filenames, block::Block)

Adds a block to the specified filenames.
"""
function setinputinfo!(job::DFJob, filenames, data::InputInfo)
    UNDO_JOBS[job.id] = deepcopy(job)

    for input in inputs(job, filenames)
        setinputinfo!(input, data)
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
            warning("Pseudo for $(id(atom)) at index $i not found in set $pseudo_set.")
        else
            pseudo!(atom, pseudo)
        end
    end
    setflags!(job, :pseudo_dir => "'$(getdefault_pseudodir(pseudo_set))'")
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

"""
    addwancalc!(job::DFJob, k_grid;
                       nscf_file          = "nscf.in",
                       wan_file           = "wan.win",
                       pw2wan_file        = "pw2wan.in",
                       wan_run_command    = "~/bin/wannier90.x ",
                       pw2wan_run_command = "mpirun -np 24 ~/bin/pw2wannier90.x",
                       wan_flags          = nothing,
                       pw2wan_flags       = nothing)

Adds a wannier calculation to a job. For now only works with QE.
"""
function addwancalc!(job::DFJob, k_grid;
                       nscf_file          = "nscf.in",
                       wan_file           = "wan.win",
                       pw2wan_file        = "pw2wan.in",
                       wan_run_command    = Exec("wannier90.x", "~/bin/"),
                       pw2wan_run_command = Exec("mpirun", "", Dict(:np => 24)),
                       pw2wan_exec        = Exec("pw2wannier90.x","~/bin/"),
                       inner_window       = (0., 0.), #no window given
                       outer_window       = (0., 0.), #no outer window given
                       wan_flags          = Dict{Symbol, Any}(),
                       pw2wan_flags       = Dict{Symbol, Any}(),
                       projections        = nothing,
                       spin               = false)

    UNDO_JOBS[job.id] = deepcopy(job)


    if inner_window != (0., 0.)
        wan_flags = merge!(wan_flags, Dict(:dis_froz_min => inner_window[1], :dis_froz_max => inner_window[2]))
    end
    if outer_window != (0., 0.)
        wan_flags = merge!(wan_flags, Dict(:dis_win_min => outer_window[1], :dis_win_max => outer_window[2]))
    end

    nscf_calc   = nothing
    scf_calc    = nothing
    pw2wan_calc = nothing
    for calc in job.calculations
        if eltype(calc) == QE
            calculation = flag(calc, :calculation)
            if calculation == "'scf'"
                scf_calc    = calc
            elseif calculation == "'nscf'"
                nscf_calc   = calc
            elseif calc.control[1].name == :inputpp
                pw2wan_calc = calc
            end
        end
    end

    if nscf_calc != nothing
        setdata!(nscf_calc, :k_points, kgrid(k_grid..., :nscf), option=:crystal)
    elseif scf_calc != nothing
        nscf_calc = deepcopy(scf_calc)
        nscf_calc.filename = nscf_file
        setflags!(nscf_calc,:calculation => "'nscf'")
        setdata!(nscf_calc, :k_points, kgrid(k_grid..., :nscf), option=:crystal)
        setflags!(nscf_calc, :noinv=>true,:nosym=>true)
        if flag(scf_calc,:noinv) != true
            setflags!(scf_calc, :noinv =>true,:nosym=>true)
            warning("Rerun scf because noinv was not set to true, and k-points will rsult in pw2wan error.")
        end
        push!(job.calculations, nscf_calc)
    else
        error("Job needs to have at least an scf calculation.")
    end
    nscf_calc.run = true
    structure = job.structure
    std_pw2wan_flags = Dict(:prefix    => flag(scf_calc, :prefix),
                           :write_amn => true,
                           :write_mmn => true,
                           :outdir    => "'./'",
                           :seedname  => "'$(splitext(wan_file)[1])'")
    if pw2wan_flags == nothing
        pw2wan_flags = std_pw2wan_flags
    else
        pw2wan_flags = merge(std_pw2wan_flags, pw2wan_flags)
    end

    if :wannier_plot in keys(wan_flags) && wan_flags[:wannier_plot]
        pw2wan_flags[:write_unk] = true
    end

    if pw2wan_calc != nothing
        setflags!(pw2wan_calc, pw2wan_flags)
    elseif eltype(scf_calc) == QE
        pw2wan_calc = DFInput{QE}(pw2wan_file, structure, pw2wan_flags, InputInfo[], pw2wan_run_command, pw2wan_exec, true)
    end

    wan_flags[:mp_grid] = typeof(k_grid) <: Array ? k_grid : convert(Array, k_grid)

    data = [InputInfo(:kpoints, :none, kgrid(k_grid..., :wan))]

    if isempty(projections)
        for atom in job.structure.atoms
            setprojections!(atom, :random)
        end
    else
        add_projections(projections, job.structure.atoms)
    end

    wan_calc1 = DFInput{Wannier90}(wan_file, structure, wan_flags, data, wan_run_command, true, true)
    if spin
        file, ext         = splitext(pw2wan_calc.filename)
        wan_file, wan_ext = splitext(wan_calc1.filename)

        pw2wan_calc_dn = deepcopy(pw2wan_calc)
        pw2wan_calc.control[:inputpp].flags[:spin_component]    = "'up'"
        pw2wan_calc.control[:inputpp].flags[:seedname]          = "'$(wan_file * "_up")'"
        pw2wan_calc_dn.control[:inputpp].flags[:spin_component] = "'down'"
        pw2wan_calc_dn.control[:inputpp].flags[:seedname]       = "'$(wan_file * "_dn")'"

        pw2wan_calc.filename    = file * "_up" * ext
        pw2wan_calc_dn.filename = file * "_dn" * ext

        wan_calc_dn1 = deepcopy(wan_calc1)
        wan_calc1.filename    = wan_file * "_up" * wan_ext
        wan_calc_dn1.filename = wan_file * "_dn" * wan_ext

        wan_calc2 = deepcopy(wan_calc1)
        push!(job.calculations, wan_calc1)
        push!(job.calculations, pw2wan_calc)
        push!(job.calculations, wan_calc2)

        wan_calc_dn2 = deepcopy(wan_calc_dn1)
        push!(job.calculations, wan_calc_dn1)
        push!(job.calculations, pw2wan_calc_dn)
        push!(job.calculations, wan_calc_dn2)

    else
        wan_calc2 = deepcopy(wan_calc1)
        push!(job.calculations, wan_calc1)
        push!(job.calculations, pw2wan_calc)
        push!(job.calculations, wan_calc2)
    end
    return job
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

path(job::DFJob, calc_filename::String) =
    joinpath(job.local_dir, input(job, calc_filename).filename)

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
sets the projections of the specified atoms inside the job structure.
"""
setprojections!(job::DFJob, projections...) = setprojections!(job.structure, projections...)
