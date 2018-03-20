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
    function DFJob(name, structure, calculations, local_dir, server, server_dir, header = get_default_job_header())
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
                    server=get_default_server(),
                    server_dir="",
                    package=:qe,
                    bin_dir="~/bin/",
                    run_command=Exec("mpirun", bin_dir, Dict(:np => 24)),
                    pseudo_set=:default,
                    pseudo_specifier="",
                    header=get_default_job_header())

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
    for (calc, data) in calculations
        calc_ = typeof(calc) == String ? Symbol(calc) : calc
        if in(calc_, [Symbol("vc-relax"), :relax, :scf])
            k_points = get(data, :k_points, [1, 1, 1, 0, 0, 0])
            k_option = :automatic
        elseif calc_ == :nscf
            k_points = get(data, :k_points, (1, 1, 1))
            k_grid   = gen_k_grid(k_points..., :nscf)
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
        push!(flags, :calculation => "'$(string(calc_))'")
        input_ = QEInput(string(calc_) * ".in",
                         QEControlBlock[],
                         [QEDataBlock(:k_points, k_option, k_points)],
                         run_command, Exec("pw.x", bin_dir),
                         true)
        set_flags!(input_, req_flags..., print=false)
        set_flags!(input_, flags..., print=false)
        push!(job_calcs, input_)
    end
    out = DFJob(job_name, structure, job_calcs, local_dir, server, server_dir, header)
    change_atoms!(out, structure.atoms, pseudo_set = pseudo_set, pseudo_specifier= pseudo_specifier)
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
        change_cell!(newjob, cell)
    end
    if atoms != nothing
        if pseudo_set == nothing
            pseudo_set, specifier = getpseudoset(job.structure.atoms[1])
            specifier = pseudo_specifier == nothing ? specifier : pseudo_specifier
            change_atoms!(newjob, atoms, pseudo_set = pseudo_set, pseudo_specifier=specifier)
        else
            change_atoms!(newjob, atoms, pseudo_set = pseudo_set, pseudo_specifier= pseudo_specifier)
        end
    end
    if pseudo_set != nothing && atoms == nothing
        change_pseudo_set!(newjob, pseudo_set, pseudo_specifier)
    end
    if server_dir != nothing
        change_server_dir!(newjob, server_dir)
    end
    if local_dir != nothing
        change_local_dir!(newjob, local_dir)
    end
    if name != nothing
        newjob.name = name
    end

    stfls!(newjob, flagstoset...)
    return newjob
end



#-------------------BEGINNING GENERAL SECTION-------------#
#all get_inputs return arrays, get_input returns the first element if multiple are found
"""
    get_inputs(job::DFJob, filenames::Array)

Returns an array of the inputs that match one of the filenames.
"""
function get_inputs(job::DFJob, filenames::Vector)
    out = DFInput[]
    for name in filenames
        push!(out, filter(x -> contains(x.filename, name), job.calculations)...)
    end
    return out
end

"""
    get_inputs(job::DFJob, filename::String)

Returns an array of the input that matches the filename.
"""
function get_inputs(job::DFJob, filename::String)
    return filter(x -> contains(x.filename, filename), job.calculations)
end

"""
    get_input(job::DFJob, filename::String)

Returns the input that matches the filename.
"""
function get_input(job::DFJob, filename::String)
    return getfirst(x -> contains(x.filename, filename), job.calculations)
end

"""
    get_input(job::DFJob, filenames::Array)

Returns an array of the inputs that match one of the filenames.
"""
function get_input(job::DFJob, filenames::Vector{String})
    return get_inputs(job, filenames)
end

"""
    load_job(job_dir::String, T=Float64; job_fuzzy = "job", new_job_name=nothing, new_homedir=nothing, server=get_default_server(),server_dir="")

Loads and returns a DFJob. If local_dir is not specified the job directory will be registered as the local one.
"""
function load_job(job_dir::String, T=Float64;
                  job_fuzzy     = "job",
                  new_job_name  = "",
                  new_local_dir = nothing,
                  server        = get_default_server(),
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

#TODO should we also create a config file for each job with stuff like server etc? and other config things,
#      which if not supplied could contain the default stuff?
"""
    pull_job(server::String, server_dir::String, local_dir::String; job_fuzzy="*job*")

Pulls job from server. If no specific inputs are supplied it pulls all .in and .tt files.
"""
function pull_job(server::String, server_dir::String, local_dir::String; job_fuzzy="*job*")
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

pull_job(args...; kwargs...) = pull_job(get_default_server(), args..., kwargs...)


"""
    load_server_job(server_dir::String, local_dir::String; server=get_default_server(), job_fuzzy="*job*", new_job_name="")

Pulls a server job to local directory and then loads it. A fuzzy search for the job file will be performed and the found input files will be pulled.
"""
function load_server_job(server_dir::String, local_dir::String;
                         server = get_default_server(),
                         job_fuzzy    = "*job*",
                         new_job_name = "")

    pull_job(server, server_dir, local_dir)
    return load_job(local_dir, server=server, server_dir=server_dir, new_job_name=new_job_name)
end

"""
    save_job(job::DFJob)

Saves a DFJob, it's job file and all it's input files.
"""
function save_job(job::DFJob)
    local_dir = job.local_dir
    if local_dir == ""
        error("Please specify a valid local_dir!")
    end
    local_dir = form_directory(job.local_dir)
    if !ispath(local_dir)
        mkpath(local_dir)
    end
    job.local_dir = local_dir
    return write_job_files(job)
end

#Incomplete everything is hardcoded for now still!!!! make it configurable
"""
    push_job(job::DFJob)

Pushes a DFJob from it's local directory to its server side directory.
"""
function push_job(job::DFJob, newfiles)
    if !ispath(job.local_dir)
        error("Please save the job locally first using save_job(job)!")
    else
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
end

#TODO only uses qsub for now. how to make it more general?
"""
    submit_job(job::DFJob; server=nothing, server_dir=nothing)

Submit a DFJob. First saves it locally, pushes it to the server then runs the job file on the server.
"""
function submit_job(job::DFJob; server=nothing, server_dir=nothing)
    if job.server == "" && server == nothing
        error("Please specify a valid server address!")
    elseif job.server_dir == "" && server_dir == nothing
        error("Please specify a valid server directory!")
    end

    if server != nothing
        job.server = server
    end
    if server_dir != nothing
        job.server_dir = server_dir
    end

    new_files = save_job(job)
    push_job(job, new_files)
    run(`ssh -t $(job.server) cd $(job.server_dir) '&&' qsub job.tt`)
end

"""
    add_calculation!(job::DFJob, input::DFInput, index::Int=length(job.calculations)+1; run_command=input.run_command, filename=input.filename)

Adds a calculation to the job, at the specified index.
"""
function add_calculation!(job::DFJob, input::DFInput, index::Int=length(job.calculations)+1;
                          run_command = input.run_command,
                          filename    = input.filename)

    UNDO_JOBS[job.id] = deepcopy(job)

    input.filename = filename
    input.run_command = run_command
    insert!(job.calculations, index, input)
    return job
end

"""
    change_flags!(job::DFJob, new_flag_data...)

Looks through all the calculations for the specified flags. If any that match and have the same types are found, they will get replaced by the new ones.
"""
function change_flags!(job::DFJob, new_flag_data...)
    UNDO_JOBS[job.id] = deepcopy(job)

    calc_filenames = [calc.filename for calc in job.calculations]
    change_flags!(job, calc_filenames, new_flag_data...)
    return job
end

"""
    change_flags!(job::DFJob, calc_filenames::Vector{String}, new_flag_data...)

Looks through the given calculations for the specified flags. If any that match and have the same types are found, they will get replaced by the new ones.
"""
function change_flags!(job::DFJob, calc_filenames::Vector{String}, new_flag_data...)
    UNDO_JOBS[job.id] = deepcopy(job)

    found_keys = Symbol[]
    for calc in get_inputs(job, calc_filenames)
        t_found_keys, = change_flags!(calc, new_flag_data...)
        for key in t_found_keys
            if !(key in found_keys) push!(found_keys, key) end
        end
    end

    n_found_keys = Symbol[]
    for (k, v) in new_flag_data
        if !(k in found_keys) push!(n_found_keys, k) end
    end

    if 1 < length(n_found_keys)
        dfprintln("flags '$(join(":" .* String.(n_found_keys),", "))' were not found in any input file, please set them first!")
    elseif length(n_found_keys) == 1
        dfprintln("flag '$(":"*String(n_found_keys[1]))' was not found in any input file, please set it first!")
    end
    return job
end
change_flags!(job::DFJob, filename::String, args...) = change_flags!(job, [filename], args...)

"""
    set_flags!(job::DFJob, calculations::Vector{String}, flags...; print=true)

Sets the flags in the calculations to the flags specified.
This only happens if the specified flags are valid for the calculations.
If necessary the correct control block will be added to the calculation (e.g. for QEInputs).

The values that are supplied will be checked whether they are valid.
"""
function set_flags!(job::DFJob, calculations::Vector{String}, flags...; print=true)
    UNDO_JOBS[job.id] = deepcopy(job)

    found_keys = Symbol[]
    for calc in get_inputs(job, calculations)
        t_found_keys, = set_flags!(calc, flags..., print=print)
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
set_flags!(job::DFJob, flags...)                   = set_flags!(job, [calc.filename for calc in job.calculations], flags...)
set_flags!(job::DFJob, filename::String, flags...) = set_flags!(job, [filename], flags...)

"""
    get_flag(job::DFJob, calc_filenames, flag_name::Symbol)

Looks through the calculation filenames and returns the value of the specified flag.
"""
function get_flag(job::DFJob, calc_filenames, flag_name::Symbol)
    for calc in get_inputs(job, calc_filenames)
        flag = get_flag(calc, flag_name)
        if flag != nothing
            return flag
        end
    end
    error("Flag $flag_name not found in any input files.")
end

"""
    get_flag(job::DFJob, flag_name::Symbol)

Looks through all the calculations and returns the value of the specified flag.
"""
function get_flag(job::DFJob, flag_name::Symbol)
    for calc in job.calculations
        flag = get_flag(calc, flag_name)
        if flag != nothing
            return flag
        end
    end
    error("Flag $flag_name not found in any input files.")
end

#TODO Change so calculations also have a name.
#TODO change after implementing k_point change so you don't need to specify all this crap
"""
    get_data(job::DFJob, calc_filenames, block_symbol::Symbol)

Looks through the calculation filenames and returns the data with the specified symbol.
"""
function get_data(job::DFJob, calc_filenames, block_symbol::Symbol)
    for calc in get_inputs(job, calc_filenames)
        return get_data(calc, block_symbol)
    end
end

"""
    get_block(job::DFJob, calc_filenames, block_symbol::Symbol)

Looks through the calculation filenames and returns the block with the specified symbol.
"""
function get_block(job::DFJob, calc_filenames, block_symbol::Symbol)
    for calc in get_inputs(job, calc_filenames)
        return get_block(calc, block_symbol)
    end
end

"""
    change_data!(job::DFJob, calc_filenames, data_block_name::Symbol, new_block_data; option=nothing)

Looks through the calculation filenames and changes the data of the datablock with `data_block_name` to `new_block_data`.
if option is specified it will set the block option to it.
"""
function change_data!(job::DFJob, calc_filenames, data_block_name::Symbol, new_block_data; option=nothing)
    UNDO_JOBS[job.id] = deepcopy(job)

    for calc in get_inputs(job, calc_filenames)
        change_data!(calc, data_block_name, new_block_data, option=option)
    end
end

"""
    remove_flags!(job::DFJob, calc_filenames, flags...)

Looks through the calculation filenames and removes the specified flags.
"""
function remove_flags!(job::DFJob, calc_filenames::Vector{<:AbstractString}, flags...)
    UNDO_JOBS[job.id] = deepcopy(job)

    for calc in get_inputs(job, calc_filenames)
        remove_flags!(calc, flags...)
    end
    return job
end
remove_flags!(job::DFJob, filename::String, flags...) = remove_flags!(job, [filename], flags...)

"""
    remove_flags!(job::DFJob, flags...)

Looks through all the calculations and removes the flags.
"""
function remove_flags!(job::DFJob, flags...)
    UNDO_JOBS[job.id] = deepcopy(job)

    for calc in job.calculations
        remove_flags!(calc, flags...)
    end
    return job
end

"""
    set_flow!(job::DFJob, should_runs::Vector{Bool})

Sets whether calculations should be ran or not. should_runs should have the same length as the amount of calculations in the job.
"""
function set_flow!(job::DFJob, should_runs::Vector{Bool})
    UNDO_JOBS[job.id] = deepcopy(job)

    assert(length(should_runs) == length(job.calculations))

    for (calc, should_run) in zip(job.calculations, should_runs)
        calc.run = should_run
    end
    return job
end

"""
    change_flow!(job::DFJob, should_runs...)

Sets whether or not calculations should be run. Calculations are specified using their indices.
"""
function change_flow!(job::DFJob, should_runs...)
    UNDO_JOBS[job.id] = deepcopy(job)

    for (filename, run) in should_runs
        for input in get_inputs(job, filename)
            input.run = run
        end
    end
    return job
end

"""
    change_flow!(job::DFJob, filenames::Array{String,1}, should_run)

Goes throug the calculation filenames and sets whether it should run or not.
"""
function change_flow!(job::DFJob, filenames::Vector{String}, should_run)
    UNDO_JOBS[job.id] = deepcopy(job)

    for calc in get_inputs(job, filenames)
        calc.run = should_run
    end
    return job
end

"""
    change_flow!(job::DFJob, should_runs::Union{Dict{String,Bool},Array{Tuple{String,Bool}}})

Runs through the calculation filenames and sets whether it should run or not.
"""
function change_flow!(job::DFJob, should_runs::Union{Dict{String,Bool}, Vector{Tuple{String,Bool}}})
    UNDO_JOBS[job.id] = deepcopy(job)

    for (name, should_run) in should_runs
        change_flow!(job, name, should_run)
    end

    return job
end

"""
    change_run_command!(job::DFJob, inputnames, run_command::Exec)

Goes through the job calculations and if it contains one of the inputnames it sets the run command of the calculation.
"""
function change_run_command!(job::DFJob, inputnames, run_command::Exec)
    UNDO_JOBS[job.id] = deepcopy(job)

    for calc in get_inputs(job, inputnames)
        calc.run_command = run_command
        dfprintln("Run command of file '$(calc.filename)' is now: '$(calc.run_command)'")
    end

    return job
end

"""
    get_run_command(job::DFJob, inputname)

Returns the `run_command`.
"""
function get_run_command(job::DFJob, inputname)
    for calc in get_inputs(job, inputname)
        return calc.run_command
    end
end

"""
    set_runflags!(job::DFJob, inputnames, flags...)

Goes through the calculations of the job and if the name contains any of the `inputnames` it sets the run_command flags to the specified ones.
"""
function set_runflags!(job::DFJob, inputnames, flags...)
    calcs = get_inputs(job, inputnames)
    for calc in calcs
        for (f,v) in flags
            calc.run_command.flags[f] = v
        end
        dfprintln("run flags of calculation $(calc.filename) are now $(calc.run_command.flags).")
    end
end

"Returns the run_command flags."
get_runflags(job::DFJob, inputname) = get_input(job, inputname).run_command[2]

"""
    set_execflags!(job::DFJob, inputnames, flags...)

Goes through the calculations of the job and if the name contains any of the `inputnames` it sets the exec flags to the specified ones.
"""
function set_execflags!(job::DFJob, inputnames, flags...)
    calcs = get_inputs(job, inputnames)
    for calc in calcs
        for (f,v) in flags
            calc.exec.flags[f] = v
        end
        dfprintln("run flags of calculation $(calc.filename) are now $(calc.exec.flags).")
    end
end

"Returns the run_command flags."
get_execflags(job::DFJob, inputname) = get_input(job, inputname).exec.flags

"""
    add_block!(job::DFJob, filenames, block::Block)

Adds a block to the specified filenames.
"""
function add_block!(job::DFJob, filenames, block::Block)
    UNDO_JOBS[job.id] = deepcopy(job)

    for input in get_inputs(job, filenames)
        add_block!(input, block)
    end
    return job
end

"""
    add_data!(job::DFJob, filenames, block_symbol, data, option=:none)

Adds a block to the specified filenames.
"""
function add_data!(job::DFJob, filenames, block_symbol, data, option=:none)
    UNDO_JOBS[job.id] = deepcopy(job)

    for input in get_inputs(job, filenames)
        add_data!(input, block_symbol, data, option)
    end
    return job
end

"""
    change_filename!(job::DFJob, old_filename::String, new_filename::String)

Changes the filename from the old to the new one.
"""
function change_filename!(job::DFJob, old_filename::String, new_filename::String)
    input          = get_input(job, old_filename)
    old_filename   = input.filename
    input.filename = new_filename
    dfprintln("Input filename in job $(job.name) has been changed from '$old_filename' to '$new_filename'.")
end

#---------------------------------END GENERAL SECTION ------------------#

"""
    change_atoms!(job::DFJob, atoms::Dict{Symbol,<:Array{<:Point3,1}}, pseudo_set_name=:default, pseudo_specifier=nothing, option=:angstrom)

Sets the data blocks with atomic positions to the new one. This is done for all calculations in the job that have that data.
If default pseudopotentials are defined, a set can be specified, together with a fuzzy that distinguishes between the possible multiple pseudo strings in the pseudo set.
These pseudospotentials are then set in all the calculations that need it.
All flags which specify the number of atoms inside the calculation also gets set to the correct value.
"""
change_atoms!(job::DFJob, atoms::Dict{Symbol,<:Vector{<:Point3}}; kwargs...) = change_atoms(job, convert_2atoms(atoms); kwargs...)

function change_atoms!(job::DFJob, atoms::Vector{<:AbstractAtom}; pseudo_set = :default, pseudo_specifier="")
    UNDO_JOBS[job.id] = deepcopy(job)

    job.structure.atoms = atoms
    change_pseudo_set!(job, pseudo_set, pseudo_specifier)
    return job
end

"""
    get_atoms(job::DFJob)

Returns a list the atoms in the structure.
"""
get_atoms(job::DFJob) = job.structure.atoms
get_cell(job::DFJob)  = job.structure.cell
#automatically sets the cell parameters for the entire job, implement others
"""
    change_cell_parameters!(job::DFJob, cell_param::Matrix)

Changes the cell parameters of the structure in the job.
"""
function change_cell!(job::DFJob, cell_param::Matrix)
    UNDO_JOBS[job.id] = deepcopy(job)

    job.structure.cell = cell_param
    return job
end

"""
    change_kpoints!(job::DFJob, calc_filename, k_points)

Changes the data in the k point `DataBlock` inside the specified calculation.
"""
function change_kpoints!(job::DFJob, calc_filenames, k_points)
    UNDO_JOBS[job.id] = deepcopy(job)
    for calc in get_inputs(job, calc_filenames)
        change_kpoints!(calc, k_points)
    end
    return job
end

"""
    change_data_option!(job::DFJob, filenames::Array{String,1}, block_symbol::Symbol, option::Symbol)

Changes the option of specified data block in the specified calculations.
"""
function change_data_option!(job::DFJob, filenames::Vector{String}, block_symbol::Symbol, option::Symbol)
    UNDO_JOBS[job.id] = deepcopy(job)

    for calc in get_inputs(job, filenames)
        change_data_option!(calc, block_symbol, option)
    end
    return job
end
change_data_option!(job::DFJob, filename::String, block_symbol::Symbol, option::Symbol) = change_data_option!(job, [filename], block_symbol, option)

"""
    change_data_option!(job::DFJob, block_symbol::Symbol, option::Symbol)

Changes the option of specified data block in all calculations that have the block.
"""
change_data_option!(job::DFJob, block_symbol::Symbol, option::Symbol) = change_data_option!(job, [i.filename for i in job.calculations], block_symbol, option)

"Changes the pseudopotentials to the specified one in the default pseudo_set."
function change_pseudo_set!(job::DFJob, pseudo_set, pseudo_specifier="")
    for (i, atom) in enumerate(job.structure.atoms)
        pseudo = get_default_pseudo(atom.id, pseudo_set, pseudo_specifier=pseudo_specifier)
        atom.pseudo = pseudo == nothing ? warning("Pseudo for $(atom.id) at index $i not found in set $pseudo_set.") : pseudo
    end
    set_flags!(job, :pseudo_dir => "'$(get_default_pseudo_dir(pseudo_set))'")
    return job
end

"""
    replace_header_word!(job::DFJob, word::String, new_word::String)


Replaces the specified word in the header with the new word.
"""
function change_header_word!(job::DFJob, word::String, new_word::String)
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
    get_errors(job::DFJob)

Prints the possible error messages in outputs of the `DFJob`.
"""
function get_errors(job::DFJob)
    outputs = pull_outputs(job)
    errors  = Dict{String, Vector{String}}()
    for out in outputs
        errors[out] = read_errors(out)
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
    add_wan_calc!(job::DFJob, k_grid;
                       nscf_file          = "nscf.in",
                       wan_file           = "wan.win",
                       pw2wan_file        = "pw2wan.in",
                       wan_run_command    = "~/bin/wannier90.x ",
                       pw2wan_run_command = "mpirun -np 24 ~/bin/pw2wannier90.x",
                       wan_flags          = nothing,
                       pw2wan_flags       = nothing)

Adds a wannier calculation to a job. For now only works with QE.
"""
function add_wan_calc!(job::DFJob, k_grid;
                       nscf_file          = "nscf.in",
                       wan_file           = "wan.win",
                       pw2wan_file        = "pw2wan.in",
                       wan_run_command    = "~/bin/wannier90.x",
                       pw2wan_run_command = "mpirun -np 24",
                       pw2wan_exec        = "~/bin/pw2wannier90.x",
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
        if typeof(calc) == QEInput
            calculation = get_flag(calc, :calculation)
            if calculation == "'scf'"
                scf_calc    = calc
            elseif calculation == "'nscf'"
                nscf_calc   = calc
            elseif calc.control_blocks[1].name == :inputpp
                pw2wan_calc = calc
            end
        end
    end

    if nscf_calc != nothing
        change_data!(nscf_calc, :k_points, gen_k_grid(k_grid[1], k_grid[2], k_grid[3], :nscf), option=:crystal)
    elseif scf_calc != nothing
        nscf_calc = deepcopy(scf_calc)
        nscf_calc.filename = nscf_file
        change_flags!(nscf_calc,:calculation => "'nscf'")
        change_data!(nscf_calc, :k_points, gen_k_grid(k_grid[1], k_grid[2], k_grid[3], :nscf), option=:crystal)
        set_flags!(nscf_calc, :noinv=>true,:nosym=>true)
        if get_flag(scf_calc,:noinv) != true
            set_flags!(scf_calc, :noinv =>true,:nosym=>true)
            warning("Rerun scf because noinv was not set to true, and k-points will rsult in pw2wan error.")
        end
        push!(job.calculations, nscf_calc)
    else
        error("Job needs to have at least an scf calculation.")
    end
    nscf_calc.run = true
    structure = job.structure
    std_pw2wan_flags = Dict(:prefix    => get_flag(scf_calc, :prefix),
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
        change_flags!(pw2wan_calc, pw2wan_flags)
    elseif typeof(scf_calc) == QEInput
        pw2wan_calc = QEInput(pw2wan_file,structure, [QEControlBlock(:inputpp, pw2wan_flags)], QEDataBlock[], pw2wan_run_command, pw2wan_exec, true)
    end

    wan_flags[:mp_grid] = typeof(k_grid) <: Array ? k_grid : convert(Array, k_grid)

    data_blocks = [WannierDataBlock(:kpoints, :none, gen_k_grid(k_grid[1], k_grid[2], k_grid[3], :wan))]

    if isempty(projections)
        for atom in job.structure.atoms
            atom.projections = :random
        end
    else
        add_projections(projections, job.structure.atoms)
    end

    wan_calc1 = WannierInput(wan_file, structure, wan_flags, data_blocks, wan_run_command, true, true)
    if spin
        file, ext         = splitext(pw2wan_calc.filename)
        wan_file, wan_ext = splitext(wan_calc1.filename)

        pw2wan_calc_dn = deepcopy(pw2wan_calc)
        pw2wan_calc.control_blocks[:inputpp].flags[:spin_component]    = "'up'"
        pw2wan_calc.control_blocks[:inputpp].flags[:seedname]          = "'$(wan_file * "_up")'"
        pw2wan_calc_dn.control_blocks[:inputpp].flags[:spin_component] = "'down'"
        pw2wan_calc_dn.control_blocks[:inputpp].flags[:seedname]       = "'$(wan_file * "_dn")'"

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

Undos the last change to the calculations of the job.
"""
function undo!(job::DFJob)
    job.calculations[:] = UNDO_JOBS[job.id].calculations[:]
    dfprintln("Restored the calculations of job '$(job.name)' to their previous state.")
    return job
end

"""
    undo(job::DFJob)

Undos the last change to the calculations of the job and returns as a new one.
"""
function undo(job::DFJob)
    return deepcopy(UNDO_JOBS[job.id])
end

"""
    add_bands_calculation!(job::DFJob, k_path::Vector{Vector{T}}; filename="bands.in", run=true) where T<:AbstractFloat

Checks if there is an scf calculation in the job and takes it's inputs to generate a bands calculation along the given k-path.
"""
function add_bands_calculation!(job::DFJob, k_path::Vector{Vector{T}}; filename="bands.in", run=true) where T<:AbstractFloat
    calc = getfirst(x -> typeof(x) == QEInput && get_flag(x, :calculation) == "'scf'", job.calculations)
    bands_calc = QEInput(calc, filename, run=run, k_points=(:crystal_b, k_path))
    push!(job.calculations, bands_calc)
    return job
end

get_path(job::DFJob, calc_filename::String) =
    joinpath(job.local_dir, get_input(job, calc_filename).filename)

"""
Sets the server dir of the job.
"""
function set_server_dir!(job, dir)
    dir = form_directory(dir)
    job.server_dir = dir
    return job
end

"""
Sets the local dir of the job.
"""
function set_local_dir!(job, dir)
    dir = form_directory(dir)
    job.local_dir = dir
    return job
end

change_server_dir!(job, dir) = set_server_dir!(job, dir)
change_local_dir!(job, dir)  = set_local_dir!(job, dir)

"""
Changes the projections of the specified atoms inside the job structure.
"""
change_projections!(job::DFJob, projections...) = change_projections!(job.structure, projections...)