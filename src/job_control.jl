#-------------------BEGINNING GENERAL SECTION-------------#
#all get_inputs return arrays, get_input returns the first element if multiple are found
"""
    get_inputs(job::DFJob, filenames::Array)

Returns an array of the inputs that match one of the filenames.
"""
function get_inputs(job::DFJob, filenames::Array)
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
    return filter(x -> contains(x.filename, filename), job.calculations)[1]
end

"""
    get_input(job::DFJob, filenames::Array)

Returns an array of the inputs that match one of the filenames.
"""
function get_input(job::DFJob, filenames::Array{String,1})
    return get_inputs(job, filenames)
end

"""
    create_job(job_name, local_dir, args...; server=get_default_server(),server_dir="")

Creates a new DFJob. 
"""
function create_job(job_name, local_dir, args...; server=get_default_server(), server_dir="")
    local_dir = form_directory(local_dir)
    inputs    = DFInput[]
    for arg in args
        push!(inputs,arg)
    end
    return DFJob(job_name, inputs, local_dir, server, server_dir)
end

"""
    load_job(job_dir::String, T=Float64; job_fuzzy = "job", new_job_name=nothing, new_homedir=nothing, server=get_default_server(),server_dir="")

Loads and returns a DFJob. If local_dir is not specified the job directory will ge registered as the local one.
"""
function load_job(job_dir::String, T=Float64;
                  job_fuzzy     = "job",
                  new_job_name  = nothing,
                  new_local_dir = nothing,
                  server        = get_default_server(),
                  server_dir    = "")
    
    job_dir = form_directory(job_dir)
    
    job_data = read_job_file(job_dir * search_dir(job_dir, job_fuzzy)[1])
    filenames    = String[]
    run_commands = String[]
    should_run   = Bool[]
    if new_job_name != nothing
        job_name = new_job_name
    elseif haskey(job_data,:name)
        job_name = job_data[:name]
    else
        job_name = "noname"
    end

    for (i, file) in enumerate(job_data[:input_files])
        if length(search_dir(job_dir, file)) == 0 && job_data[:should_run][i]
            error("Error: there are calculations that should run but have no input file ($file).")
        elseif length(search_dir(job_dir, file)) != 0
            push!(filenames,    file)
            push!(run_commands, job_data[:run_commands][i])
            push!(should_run,   job_data[:should_run][i])
        end
    end

    t_calcs = Array{DFInput,1}()
    for (filename, run_command, run) in zip(filenames, run_commands, should_run)
        filename = job_dir * filename
        if contains(run_command, "wan") && !contains(run_command, "pw2wannier90")
            s_run_command = split(run_command)
            if "-pp" in s_run_command
                run_command = join(s_run_command[1:end-1])
                push!(t_calcs, read_wannier_input(filename * ".win", T,
                                                  run_command = run_command,
                                                  run         = run,
                                                  preprocess  = true))
            else
                push!(t_calcs,read_wannier_input(filename * ".win", T, 
                                                 run_command = run_command,
                                                 run         = run,
                                                 preprocess  = false))
            end
        elseif contains(run_command, "abinit")
            push!(t_calcs, read_abi_input(filename, T,
                                          run_command = run_command,
                                          pseudos     = job_data[:abinit_pseudos])...)
        else
            _t = split(run_command)
            push!(t_calcs, read_qe_input(filename, T, run_command=join(_t[1:end-1]," "), exec=_t[end], run=run))
        end
    end

    if new_local_dir != nothing
        return DFJob(job_name, t_calcs, new_local_dir, server, server_dir, job_data[:header])
    else
        return DFJob(job_name, t_calcs, job_dir, server, server_dir, job_data[:header])
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
        job_data = read_job_file(local_dir * job_file)
        for (file, run_command) in zip(job_data[:input_files], job_data[:run_commands])
            if !contains(file, ".") && contains(run_command, "wannier90.x")
                pull_server_file(file * ".win")
            else
                pull_server_file(file)
            end
        end
    end

end

pull_job(args...; kwargs...) = pull_job(get_default_server(), args..., kwargs...)


"""
    load_server_job(server::String, server_dir::String, local_dir::String; job_fuzzy="*job*", job_name=nothing)

Pulls a server job to local directory and then loads it. A fuzzy search for the job file will be performed and the found input files will be pulled.
"""
function load_server_job(server::String, server_dir::String, local_dir::String;
                         job_fuzzy    = "*job*",
                         new_job_name = nothing)

    pull_job(server,server_dir,local_dir)
    return load_job(local_dir, server=server, server_dir=server_dir, new_job_name=new_job_name)
end

load_server_job(args...; kwargs...) = load_server_job(get_default_server(), args..., kwargs...)

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
        for file in newfiles
            run(`scp $(job.local_dir * file) $(job.server * ":" * job.server_dir)`)
        end
        # for calc in job.calculations
        #     if calc.filename in newfiles
        #         run(`scp $(job.local_dir * calc.filename) $(job.server * ":" * job.server_dir)`)
        #     end
        # end
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
    print_info(input)
    print_flow(job)
end

"""
    change_flags!(job::DFJob, new_flag_data::OrderedDict{Symbol,<:Any})

Looks through all the calculations for the specified flags. If any that match and have the same types are found, they will get replaced by the new ones.
"""
function change_flags!(job::DFJob, new_flag_data...)
    UNDO_JOBS[job.id] = deepcopy(job)

    calc_filenames = [calc.filename for calc in job.calculations]
    change_flags!(job, calc_filenames, new_flag_data...)
end

"""
    change_flags!(job::DFJob, calc_filenames, new_flag_data::OrderedDict{Symbol,<:Any})

Looks through the given calculations for the specified flags. If any that match and have the same types are found, they will get replaced by the new ones.
"""
function change_flags!(job::DFJob, calc_filenames::Array{String,1}, new_flag_data...)
    UNDO_JOBS[job.id] = deepcopy(job)

    found_keys = Symbol[]
    for calc in get_inputs(job, calc_filenames)
        t_found_keys = change_flags!(calc, new_flag_data...)
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

end
change_flags!(job::DFJob, filename::String, args...) = change_flags!(job, [filename], args...)

"""
    set_flags!(job::DFJob, calculations::Array{String,1}, flags...)

Sets the flags in the calculations to the flags specified. 
This only happens if the specified flags are valid for the calculations.
If necessary the correct control block will be added to the calculation (e.g. for QEInputs).

The values that are supplied will be checked whether they are valid.
"""
function set_flags!(job::DFJob, calculations::Array{String,1}, flags...)
    UNDO_JOBS[job.id] = deepcopy(job)

    found_keys = Symbol[]
    for calc in get_inputs(job, calculations)
        t_found_keys = set_flags!(calc, flags...)
        for key in t_found_keys
            if !(key in found_keys) push!(found_keys, key) end
        end
    end

    n_found_keys = Symbol[]
    for (k, v) in flags
        if !(k in found_keys) push!(n_found_keys, k) end
    end

    if 1 < length(n_found_keys)
        dfprintln("flags '$(join(":" .* String.(n_found_keys),", "))' were not found in the allowed input variables of the specified calculations!")
    elseif length(n_found_keys) == 1
        dfprintln("flag '$(":" * String(n_found_keys[1]))' was not found in the allowed input variables of the specified calculations!")
    end
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
    error("Flag $flag not found in any input files.")
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
    error("Flag $flag not found in any input files.")
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
    change_data!(job::DFJob, calc_filenames, data_block_name::Symbol, new_block_data)

Looks through the calculation filenames and changes the data of the datablock with `data_block_name` to `new_block_data`
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
function remove_flags!(job::DFJob, calc_filenames::Array{<:String,1}, flags...)
    UNDO_JOBS[job.id] = deepcopy(job)

    for calc in get_inputs(job, calc_filenames)
        remove_flags!(calc, flags...)
    end
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
end

"""
    set_flow!(job::DFJob, should_runs::Array{Bool,1})

Sets whether calculations should be ran or not. should_runs should have the same length as the amount of calculations in the job.
"""
function set_flow!(job::DFJob, should_runs::Array{Bool,1})
    UNDO_JOBS[job.id] = deepcopy(job)

    assert(length(should_runs) == length(job.calculations))

    for (calc, should_run) in zip(job.calculations, should_runs)
        calc.run = should_run
    end
    print_flow(job)
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

    print_flow(job)
end

"""
    change_flow!(job::DFJob, filenames, should_run)

Goes throug the calculation filenames and sets whether it should run or not.
"""
function change_flow!(job::DFJob, filenames::Array{String,1}, should_run)
    UNDO_JOBS[job.id] = deepcopy(job)

    for calc in get_inputs(job, filenames)
        calc.run = should_run
    end
end

"""
    change_flow!(job::DFJob, should_runs::Union{OrderedDict{String,Bool},Array{Tuple{String,Bool}}})

Runs through the calculation filenames and sets whether it should run or not.
"""
function change_flow!(job::DFJob, should_runs::Union{OrderedDict{String,Bool},Array{Tuple{String,Bool}}})
    UNDO_JOBS[job.id] = deepcopy(job)

    for (name, should_run) in should_runs
        change_flow!(job, name, should_run)
    end

    print_flow(job)
end

"""
    change_run_command!(job::DFJob, filenames, run_command)

Goes through the calculation filenames and sets the run command of the calculation.
"""
function change_run_command!(job::DFJob, filenames, run_command)
    UNDO_JOBS[job.id] = deepcopy(job)

    for calc in get_inputs(job, filenames)
        calc.run_command = run_command
        dfprintln("Run command of file '$(calc.filename)' is now: '$(calc.run_command)'")
    end

end

"""
    get_run_command(job::DFJob, filename)

Returns the run command for the specified calculation.
"""
function get_run_command(job::DFJob, filename)
    for calc in get_inputs(job, filename)
        return calc.run_command
    end
end

"""
    print_run_command(job::DFJob, filenames)

Prints the run command of the specified calculations.
"""
function print_run_command(job::DFJob, filenames)
    for calc in get_inputs(job, filenames)
        dfprintln("Run command of file '$(calc.filename)' is: '$(calc.run_command)'.")
        dfprintln("")
    end
end

"""
    print_flow(job::DFJob)

Prints the calculation sequence of the job.
"""
function print_flow(job::DFJob)
    for (i, calc) in enumerate(job.calculations)
        dfprintln("$i: $(calc.filename) => runs: $(calc.run)")
    end
end

"""
    print_block(job::DFJob, block_name::Symbol)

Prints information of the specified block name of all the calculations in the job.
"""
function print_block(job::DFJob, block_name::Symbol)
    for calc in job.calculations
        if print_block(calc, block_name) dfprintln("") end
    end
end

"""
    print_block(job::DFJob, filenames, block_symbol::Symbol)

Prints the information of the block in a selected file of the job.
"""
function print_block(job::DFJob, filenames, block_symbol::Symbol)
    for calc in get_inputs(job, filenames)
        print_block(calc, block_symbol)
    end
end

"""
    print_blocks(job::DFJob, calc_filenames)

Prints information on all the blocks in the specified calculations.
"""
function print_blocks(job::DFJob, calc_filenames)
    for calc in get_inputs(job, calc_filenames)
        print_blocks(calc)
    end
end

"""
    print_blocks(job::DFJob)

Prints information of all the blocks of all the calculations in the job.
"""
function print_blocks(job::DFJob)
    for calc in job.calculations
        print_blocks(calc)
        dfprintln("#------------------------------------#")
    end
end
print_data(job::DFJob)                     = print_blocks(job)
print_data(job::DFJob, calc_filenames)     = print_blocks(job, calc_filenames)
print_data(job::DFJob, block_name::Symbol) = print_block(job, block_name)

"""
    print_info(job::DFJob, filenames::Array{String,1})

Prints general info of the job, and the specified filenames.
"""
function print_info(job::DFJob, filenames::Array{String,1})
    s = """--------------------
    DFJob:      $(job.name)
    Local_dir:  $(job.local_dir)
    Server:     $(job.server)
    Server_dir: $(job.server_dir)
    $(length(job.calculations)) calculations
    --------------------
    """
    dfprintln(s)

    for calc in get_inputs(job, filenames)
        print_info(calc)
        dfprintln(" ")
    end

end
print_info(job::DFJob)                   = print_info(job, [calc.filename for calc in job.calculations])
print_info(job::DFJob, filename::String) = print_info(job, [filename])

"""
    print_flags(job::DFJob)

Prints flags of all the calculations in the job.
"""
function print_flags(job::DFJob)
    for calc in job.calculations
        print_flags(calc)
    end
end

"""
    print_flags(job::DFJob, calc_filename::String)

Prints flags of the specified calculation.
"""
function print_flags(job::DFJob, calc_filename::String)
    for calc in get_inputs(job, calc_filename)
        print_flags(calc)
    end
end

"""
    print_flags(job::DFJob, calc_filenames::Array{String,1})

Prints the flags of the specified calculations.
"""
function print_flags(job::DFJob, calc_filenames::Array{String,1})
    for file in calc_filenames
        print_flags(job, file)
    end
end

"""
    print_flags(job::DFJob, flags)

Prints the specified flags running through all the claculations in the job.
"""
function print_flags(job::DFJob, flags::Array{Symbol,1})
    for flag in flags
        print_flag(job, flag)
    end
end

"""
    print_flag(job::DFJob, flag::Symbol)

Prints the specified flag running through all the calculations in the job.
"""
function print_flag(job::DFJob, flag::Symbol)
    for calc in job.calculations
        print_flag(calc, flag)
    end
end

"""
    add_block!(job::DFJob, filenames, block::Block)

Adds a block to the specified filenames.
"""
function add_block!(job::DFJob, filenames, block::Block)
    UNDO_JOBS[job.id] = deepcopy(job)

    for input in get_inputs(job, filenames)
        add_block!(input, block)
    end
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
    change_atoms!(job::DFJob, atoms::OrderedDict{Symbol,<:Array{<:Point3D,1}}, pseudo_set_name=:default, pseudo_fuzzy=nothing, option=:angstrom)

Sets the data blocks with atomic positions to the new one. This is done for all calculations in the job that have that data. 
If default pseudopotentials are defined, a set can be specified, together with a fuzzy that distinguishes between the possible multiple pseudo strings in the pseudo set. 
These pseudospotentials are then set in all the calculations that need it.
All flags which specify the number of atoms inside the calculation also gets set to the correct value.
"""
function change_atoms!(job::DFJob, atoms::OrderedDict{Symbol,<:Array{<:Point3D,1}}; kwargs...)
    UNDO_JOBS[job.id] = deepcopy(job)

    for calc in job.calculations
        change_atoms!(calc, atoms; kwargs...)
    end
    
end

"""
    get_atoms(job::DFJob, calc_filename)

Returns a list of the atomic positions in Angstrom.
"""
function get_atoms(job::DFJob, calc_filename)
    return get_atoms(get_input(job, calc_filename))
end

"""
    sync_atoms!(job::DFJob, template_filename::String; kwargs...)

Syncs the atoms to the atoms specified in the template input.
"""
function sync_atoms!(job::DFJob, template_filename::String; kwargs...)
    input = get_input(job, template_filename)
    atoms = get_atoms(input)
    change_atoms!(job, atoms; kwargs...)
end

"""
    sync_cell!(job::DFJob, template_filename::String)

Syncs the cells of all input files to the one given in the template file.
All cells will be put in Angstrom coordinates.
"""
function sync_cell!(job::DFJob, template_filename::String)
    calc = get_input(job, template_filename)
    cell = get_cell(calc)
    for input in job.calculations
        change_cell!(input, cell, option=:angstrom)
    end
end

#automatically sets the cell parameters for the entire job, implement others
"""
    change_cell_parameters!(job::DFJob, cell_param::Matrix)

Changes the cell parameters in all the input files that have that `DataBlock`.
"""
function change_cell!(job::DFJob, cell_param::Matrix)
    UNDO_JOBS[job.id] = deepcopy(job)

    for calc in job.calculations
        change_cell!(calc, cell_param)
    end
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
end

"""
    change_data_option!(job::DFJob, filenames::Array{String,1}, block_symbol::Symbol, option::Symbol)

Changes the option of specified data block in the specified calculations.
"""
function change_data_option!(job::DFJob, filenames::Array{String,1}, block_symbol::Symbol, option::Symbol)
    UNDO_JOBS[job.id] = deepcopy(job)

    for calc in get_inputs(job, filenames)
        change_data_option!(calc, block_symbol, option)
    end
end
change_data_option!(job::DFJob, filename::String, block_symbol::Symbol, option::Symbol) = change_data_option!(job, [filename], block_symbol, option)

"""
    change_data_option!(job::DFJob, block_symbol::Symbol, option::Symbol)

Changes the option of specified data block in all calculations that have the block.
"""
change_data_option!(job::DFJob, block_symbol::Symbol, option::Symbol) = change_data_option!(job, [i.filename for i in job.calculations], block_symbol, option)

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

end

"""
    get_errors(job::DFJob)

Prints the possible error messages in outputs of the `DFJob`.
"""
function get_errors(job::DFJob)
    outputs = pull_outputs(job)
    errors  = OrderedDict{String,Array{String,1}}()
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
                       wan_run_command    = "~/bin/wannier90.x ",
                       pw2wan_run_command = "mpirun -np 24 ~/bin/pw2wannier90.x",
                       inner_window       = (0., 0.), #no window given 
                       outer_window       = (0., 0.), #no outer window given
                       wan_flags          = OrderedDict{Symbol, Any}(),
                       pw2wan_flags       = OrderedDict{Symbol, Any}(),
                       projections        = nothing,
                       spin               = false)
            
    UNDO_JOBS[job.id] = deepcopy(job)


    if inner_window != (0., 0.) #scalarize
        wan_flags = merge!(wan_flags, OrderedDict(:dis_froz_min => inner_window[1], :dis_froz_max => inner_window[2]))
    end
    if outer_window != (0., 0.) 
        wan_flags = merge!(wan_flags, OrderedDict(:dis_win_min => outer_window[1], :dis_win_max => outer_window[2]))
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
    
    std_pw2wan_flags = OrderedDict(:prefix    => get_flag(scf_calc, :prefix),
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
        pw2wan_calc = QEInput(pw2wan_file, [QEControlBlock(:inputpp, pw2wan_flags)], QEDataBlock[], pw2wan_run_command, true)
    end

    cell_block = get_block(scf_calc, :cell_parameters)
    if cell_block.option == :alat
        alat = 1.0
        if get_flag(scf_calc, :A) != nothing
            alat = get_flag(scf_calc, :A)
        elseif get_flag(scf_calc, Symbol("celldm(1)")) != nothing
            alat = conversions[:bohr2ang] * get_flag(scf_calc, Symbol("celldm(1)"))
        else
            error("Please set either flag :A or :celldm(1) when cell_parameters are in alat.")
        end
        cell = cell_block.data .* alat
    elseif cell_block.option == :bohr
        cell = cell_block.data .* conversions[:bohr2ang]
    else
        cell = cell_block.data
    end
    
    atoms_block = get_block(scf_calc, :atomic_positions)
    atoms = deepcopy(atoms_block.data)
    if atoms_block.option == :alat
        alat = 1.0
        if get_flag(scf_calc, :A) != nothing
            alat = get_flag(scf_calc, :A)
        elseif get_flag(scf_calc, Symbol("celldm(1)")) != nothing
            alat = conversions[:bohr2ang] * get_flag(scf_calc, Symbol("celldm(1)"))
        else
            error("Please set either flag :A or :celldm(1) when cell_parameters are in alat.")
        end
        atoms[key] = positions .* alat
    elseif atoms_block.option == :bohr
        for (key, positions) in atoms_block.data
            atoms[key] = positions .* conversions[:bohr2ang]
        end
    elseif atoms_block.option == :crystal
        for (key, positions) in atoms_block.data
            t_pos = Point3D[]
            for p in positions
                push!(t_pos, cell * p)
            end
            atoms[key] = t_pos
        end
    end
    wan_flags[:mp_grid] = typeof(k_grid) <: Array ? k_grid : convert(Array, k_grid)

    data_blocks = [WannierDataBlock(:atoms_cart, :ang, atoms),
                   WannierDataBlock(:unit_cell_cart, :ang, cell), 
                   WannierDataBlock(:kpoints, :none, gen_k_grid(k_grid[1], k_grid[2], k_grid[3], :wan))]

    if isempty(projections)
        push!(data_blocks, WannierDataBlock(:projections, :random, nothing))
    else
        push!(data_blocks, WannierDataBlock(:projections, :none, projections))
    end

    wan_calc1 = WannierInput(wan_file, wan_flags, data_blocks, wan_run_command, true, true)
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
    print_info(job)
end

"""
    undo!(job::DFJob)

Undos the last change to the calculations of the job.
"""
function undo!(job::DFJob)
    job.calculations[:] = UNDO_JOBS[job.id].calculations[:]
    dfprintln("Restored the calculations of job '$(job.name)' to their previous state.")
end

"""
    undo(job::DFJob)

Undos the last change to the calculations of the job and returns as a new one.
"""
function undo(job::DFJob)
    return deepcopy(UNDO_JOBS[job.id])
end

"""
    add_bands_calculation!(job::DFJob, k_path::Array{Array{<:AbstractFloat,1},1})

Checks if there is an scf calculation in the job and takes it's inputs to generate a bands calculation along the given k-path.
"""
function add_bands_calculation!(job::DFJob, k_path::Array{Array{T,1},1}; filename="bands.in", run=true) where T<:AbstractFloat
    for calc in job.calculations
        if typeof(calc) == QEInput && get_flag(calc,:calculation) == "'scf'"
           bands_calc = QEInput(calc, filename, run=run, k_points=(:crystal_b, k_path))
           push!(job.calculations, bands_calc)
           break
        end
    end
end

