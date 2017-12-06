
#---------------------------------BEGINNING GENERAL SECTION ---------------------#

"""
create_job(job_name, local_dir, args...; server=get_default_server(),server_dir="")

Creates a new DFJob. 
"""
function create_job(job_name, local_dir, args...; server=get_default_server(),server_dir="")
  local_dir = form_directory(local_dir)
  inputs    = DFInput[]
  for arg in args
    push!(inputs,arg)
  end
  return DFJob(job_name,inputs,local_dir,server,server_dir)
end

"""
load_job(job_dir::String, T=Float32; job_fuzzy = "job", new_job_name=nothing, new_homedir=nothing, server=get_default_server(),server_dir="")

Loads and returns a DFJob. If local_dir is not specified the job directory will ge registered as the local one.
"""
#should we call this load local job?
function load_job(job_dir::String, T=Float32; job_fuzzy = "job", new_job_name=nothing, new_local_dir=nothing, server=get_default_server(),server_dir="")
  job_dir = form_directory(job_dir)
  
  job_name,t_inputs,t_outputs,t_run_commands,t_should_run,header = read_job_file(job_dir*search_dir(job_dir,job_fuzzy)[1])
  filenames    = String[]
  run_commands = String[]
  should_run   = Bool[]
  
  if new_job_name != nothing
    job_name = new_job_name
  end
  for (i,file) in enumerate(t_inputs)
    if length(search_dir(job_dir,file))==0 && t_should_run[i]
      error("Error: there are calculations that should run but have no input file ($file).")
    elseif length(search_dir(job_dir,file))!=0
      push!(filenames,file)
      push!(run_commands,t_run_commands[i])
      push!(should_run,t_should_run[i])
    end
  end
  
  t_calcs = Array{DFInput,1}()
  for (filename,run_command,run) in zip(filenames,run_commands,should_run)
    filename = job_dir*filename
    if contains(run_command,"wan") && !contains(run_command,"pw2wannier90")
      s_run_command = split(run_command)
      if "-pp" in s_run_command
        run_command = join(s_run_command[1:end-1])
        push!(t_calcs,read_wannier_input(filename*".win",T,run_command=run_command,run=run,preprocess=true))
      else
        push!(t_calcs,read_wannier_input(filename*".win",T,run_command=run_command,run=run,preprocess=false))
      end
    else
      push!(t_calcs,read_qe_input(filename,T,run_command=run_command,run=run))
    end
  end
  if new_local_dir != nothing
    return DFJob(job_name,t_calcs,new_local_dir,server,server_dir,header)
  else
    return DFJob(job_name,t_calcs,job_dir,server,server_dir,header)
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
  pull_server_file(filename) = pull_file(server,server_dir,local_dir,filename)
  pull_server_file(job_fuzzy)
  job_file = search_dir(local_dir,strip(job_fuzzy,'*'))[1]
  if job_file != nothing
    job_name,inputs,outputs, run_commands , _ = read_job_file(local_dir*job_file)
    for (file,run_command) in zip(inputs,run_commands)
      if !contains(file,".") && contains(run_command,"wannier90.x")
        pull_server_file(file*".win")
      else
        pull_server_file(file)
      end
    end
  end
end

pull_job(args...;kwargs...) = pull_job(get_default_server(),args...,kwargs...)


"""
    load_server_job(server::String, server_dir::String, local_dir::String; job_fuzzy="*job*", job_name=nothing)

Pulls a server job to local directory and then loads it. A fuzzy search for the job file will be performed and the found input files will be pulled.
"""
function load_server_job(server::String, server_dir::String, local_dir::String; job_fuzzy="*job*", new_job_name=nothing)
  pull_job(server,server_dir,local_dir)
  return load_job(local_dir,server=server,server_dir=server_dir,new_job_name = new_job_name)
end

load_server_job(args...;kwargs...) = load_server_job(get_default_server(),args...,kwargs...)

"""
    save_job(df_job::DFJob)

Saves a DFJob, it's job file and all it's input files.
"""
function save_job(df_job::DFJob)
  local_dir = df_job.local_dir
  if local_dir == ""
    error("Please specify a valid local_dir!")
  end
  local_dir = form_directory(df_job.local_dir)
  if !ispath(local_dir)
    mkpath(local_dir)
  end
  df_job.local_dir = local_dir
  write_job_files(df_job)
end

#Incomplete everything is hardcoded for now still!!!! make it configurable
"""
    push_job(df_job::DFJob)

Pushes a DFJob from it's local directory to its server side directory.
"""
function push_job(df_job::DFJob)
  if !ispath(df_job.local_dir)
    error("Please save the job locally first using save_job(job)!")
  else
    
    for calc in df_job.calculations
      run(`scp $(df_job.local_dir*calc.filename) $(df_job.server*":"*df_job.server_dir)`)
    end
    run(`scp $(df_job.local_dir*"job.tt") $(df_job.server*":"*df_job.server_dir)`)
  end
end

#TODO only uses qsub for now. how to make it more general?
"""
    submit_job(df_job::DFJob)

Submit a DFJob. First saves it locally, pushes it to the server then runs the job file on the server.
"""
function submit_job(df_job::DFJob; server=nothing, server_dir=nothing)
  if df_job.server == "" && server == nothing
    error("Please specify a valid server address!")
  elseif df_job.server_dir == "" && server_dir == nothing
    error("Please specify a valid server directory!")
  end
  if server != nothing
    df_job.server = server
  end
  if server_dir != nothing
    df_job.server_dir = server_dir
  end
  save_job(df_job)
  push_job(df_job)
  run(`ssh -t $(df_job.server) cd $(df_job.server_dir) '&&' qsub job.tt`)
end

"""
    add_calculation!(df_job::DFJob, input::DFInput, run_index::Int=length(df_job.calculations)+1; run_command=input.run_command, filename=input.filename)

Adds a calculation to the job, at the specified run_index.
"""
function add_calculation!(df_job::DFJob, input::DFInput, run_index::Int=length(df_job.calculations)+1; run_command=input.run_command, filename=input.filename)
  input.filename = filename
  input.run_command = run_command
  insert!(df_job.calculations,run_index,input)
  print_info(input)
  print_flow(df_job)
end

"""
    change_flags!(df_job::DFJob, new_flag_data::Dict{Symbol,<:Any})

Looks through all the calculations for the specified flags. If any that match and have the same types are found, they will get replaced by the new ones.
"""
function change_flags!(df_job::DFJob, new_flag_data...)
  calc_filenames = [calc.filename for calc in df_job.calculations]
  change_flags!(df_job,calc_filenames,new_flag_data...)
end

"""
    change_flags!(df_job::DFJob, calc_filenames, new_flag_data::Dict{Symbol,<:Any})

Looks through the given calculations for the specified flags. If any that match and have the same types are found, they will get replaced by the new ones.
"""
function change_flags!(df_job::DFJob, calc_filenames::Array{String,1}, new_flag_data...)
  found_keys = Symbol[]
  for calc in get_inputs(df_job,calc_filenames)
    t_found_keys = change_flags!(calc,new_flag_data...)
    for key in t_found_keys
      if !(key in found_keys) push!(found_keys,key) end
    end
  end
  n_found_keys = Symbol[]
  for (k,v) in new_flag_data
    if !(k in found_keys) push!(n_found_keys,k) end
  end
  if 1 < length(n_found_keys)
    println("flags '$(join(":" .* String.(n_found_keys),", "))' were not found in any input file, please set them first!")
  elseif length(n_found_keys) == 1
    println("flag '$(":"*String(n_found_keys[1]))' was not found in any input file, please set it first!")
  end
end
change_flags!(df_job::DFJob, filename::String, args...) = change_flags!(df_job,[filename],args...)
"""
    get_flag(df_job::DFJob, calc_filenames, flag::Symbol)

Looks through the calculation filenames and returns the value of the specified flag.
"""
function get_flag(df_job::DFJob, calc_filenames, flag::Symbol)
  for calc in get_inputs(df_job,calc_filenames)
    return get_flag(calc,flag)
  end
end

"""
    get_flag(df_job::DFJob, flag::Symbol)

Looks through all the calculations and returns the value of the specified flag.
"""
function get_flag(df_job::DFJob, flag::Symbol)
  for calc in df_job.calculations
    return get_flag(calc,flag)
  end
end

#TODO Change so calculations also have a name.
#TODO change after implementing k_point change so you don't need to specify all this crap
"""
    get_data(df_job::DFJob, calc_filenames, block_symbol::Symbol)

Looks through the calculation filenames and returns the data with the specified symbol.
"""
function get_data(df_job::DFJob, calc_filenames, block_symbol::Symbol)
  for calc in get_inputs(df_job,calc_filenames)
    return get_data(calc,block_symbol)
  end
end

"""
    get_block(df_job::DFJob, calc_filenames, block_symbol::Symbol)

Looks through the calculation filenames and returns the block with the specified symbol.
"""
function get_block(df_job::DFJob, calc_filenames, block_symbol::Symbol)
  for calc in get_inputs(df_job,calc_filenames)
    return get_block(calc,block_symbol)
  end
end

"""
    change_data!(df_job::DFJob, calc_filenames, data_block_name::Symbol, new_block_data)

Looks through the calculation filenames and changes the data of the datablock with `data_block_name` to `new_block_data`
"""
function change_data!(df_job::DFJob, calc_filenames, data_block_name::Symbol, new_block_data)
  for calc in get_inputs(df_job,calc_filenames)
    change_data!(calc,data_block_name,new_block_data)
  end
end

#I'm not sure if this is a good idea. Maybe require explicitely also the input filename
#TODO after making defaults this can be done better
"""
   add_flags!(df_job::DFJob, calc_filenames::Array{String,1}, control_block_name::Symbol, flags)

Adds the flags to the controlblocks of the specified inputs with the control block. This assumes that there are `ControlBlocks` in the calculations e.g. in `QEInput`.
"""
function add_flags!(job::DFJob, calc_filenames::Array{String,1}, control_block_name::Symbol, flags...)
  for calc in get_inputs(job,calc_filenames)
    if :control_blocks in fieldnames(calc)
      add_flags!(calc,control_block_name,flags...)
    end
  end
end
"""
   add_flags!(df_job::DFJob, control_block_name::Symbol, flags)

Adds the flags to the controlblocks of all inputs with the control block. This assumes that there are `ControlBlocks` in the calculations e.g. in `QEInput`.
"""
add_flags!(df_job::DFJob, control_block_name::Symbol, flags...) = add_flags!(df_job,[i.filename for i in df_job.calculations], control_block_name,flags...)

add_flags!(job::DFJob,filename::String,control_block_name::Symbol,flags...)=add_flags!(job,[filename],control_block_name,flags...)

"""
    add_flags!(df_job::DFJob,filenames::Array{String,1}, flags...)

Adds the flags to specified input files. This assumes that the input has a field `flags`. Works for e.g. `WannierInput`.
"""

function add_flags!(df_job::DFJob,filenames::Array{String,1}, flags...)
  for calc in get_inputs(df_job,filenames)
    if :flags in fieldnames(calc)
      add_flags!(calc,flags...)
    end
  end
end
"""
    add_flags!(job::DFJob,flags...)

Adds the flags to all inputs that have fieldname 'flags'. Works for e.g. 'WannierInput'.
"""
add_flags!(job::DFJob,flags...) = add_flags!(job,[i.filename for i in job.calculations],flags...)
add_flags!(df_job::DFJob, filename::String,flags...) = add_flags!(df_job,[filename],flags...)
"""
    remove_flags!(df_job::DFJob, calc_filenames, flags...)

Looks through the calculation filenames and removes the specified flags.
"""
function remove_flags!(df_job::DFJob, calc_filenames::Array{<:String,1}, flags...)
  for calc in get_inputs(df_job,calc_filenames)
    remove_flags!(calc,flags...)
  end
end

remove_flags!(job::DFJob,filename::String,flags...)=remove_flags!(job,[filename],flags...)
"""
    remove_flags!(df_job::DFJob, flags...)

Looks through all the calculations and removes the flags.
"""
function remove_flags!(df_job::DFJob, flags...)
  for calc in df_job.calculations
    remove_flags!(calc,flags...)
  end
end

"""
    set_flow!(df_job::DFJob, should_runs::Array{Bool,1})

Sets whether calculations should be ran or not. should_runs should have the same length as the amount of calculations in the job.
"""
function set_flow!(df_job::DFJob, should_runs::Array{Bool,1})
  assert(length(should_runs)==length(df_job.calculations))
  for (calc,should_run) in zip(df_job.calculations,should_runs)
    calc.run = should_run
  end
  print_flow(df_job)
end

"""
    change_flow!(df_job::DFJob, should_runs...)

Sets whether or not calculations should be run. Calculations are specified using their indices.
"""
function change_flow!(df_job::DFJob, should_runs...)
  for (filename,run) in should_runs
    get_input(df_job,filename).run = run
  end
  print_flow(df_job)
end

"""
    change_flow!(df_job::DFJob, filenames, should_run)

Goes throug the calculation filenames and sets whether it should run or not.
"""
function change_flow!(df_job::DFJob, filenames, should_run)
  for calc in get_inputs(df_job,filenames)
    calc.run = should_run
  end
end

"""
    change_flow!(df_job::DFJob, should_runs::Union{Dict{String,Bool},Array{Tuple{String,Bool}}})

Runs through the calculation filenames and sets whether it should run or not.
"""
function change_flow!(df_job::DFJob, should_runs::Union{Dict{String,Bool},Array{Tuple{String,Bool}}})
  for (name,should_run) in should_runs
    change_flow!(df_job,name,should_run)
  end
  print_flow(df_job)
end

"""
    change_run_command!(df_job::DFJob, filenames, run_command)

Goes through the calculation filenames and sets the run command of the calculation.
"""
function change_run_command!(df_job::DFJob, filenames, run_command)
  for calc in get_inputs(df_job,filenames)
    calc.run_command = run_command
    println("Run command of file '$(calc.filename)' is now: '$(calc.run_command)'")
  end
end

"""
    get_run_command(df_job::DFJob, filename)

Returns the run command for the specified calculation.
"""
function get_run_command(df_job::DFJob, filename)
  for calc in get_inputs(df_job,filename)
    return calc.run_command
  end
end

"""
    print_run_command(df_job::DFJob, filenames)

Prints the run command of the specified calculations.
"""
function print_run_command(df_job::DFJob, filenames)
  for calc in get_inputs(df_job,filenames)
    println("Run command of file '$(calc.filename)' is: '$(calc.run_command)'.")
    println("")
  end
end

"""
    print_flow(df_job::DFJob)

Prints the calculation sequence of the job.
"""
function print_flow(df_job::DFJob)
  for (i,calc) in enumerate(df_job.calculations)
    println("$i: $(calc.filename) => runs: $(calc.run)")
  end
end

"""
    print_block(job::DFJob, block_name::Symbol)

Prints information of the specified block name of all the calculations in the job.
"""
function print_block(job::DFJob, block_name::Symbol)
  for calc in job.calculations
    if print_block(calc,block_name) println("") end
  end
end

"""
    print_block(job::DFJob, filenames, block_symbol::Symbol)

Prints the information of the block in a selected file of the job.
"""
function print_block(job::DFJob, filenames, block_symbol::Symbol)
  for calc in get_inputs(job,filenames)
    print_block(calc,block_symbol)
  end
end

"""
    print_blocks(job::DFJob, calc_filenames)

Prints information on all the blocks in the specified calculations.
"""
function print_blocks(job::DFJob, calc_filenames)
  for calc in get_inputs(job,calc_filenames)
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
    println("#------------------------------------#")
  end
end

print_data(job::DFJob) = print_blocks(job)
print_data(job::DFJob, calc_filenames) = print_blocks(job, calc_filenames)
print_data(job::DFJob, block_name::Symbol) = print_block(job, block_name)

"""
    print_info(job::DFJob)

Prints general info of the job.
"""
function print_info(job::DFJob)
  println("--------------------")
  println("DFJob:      $(job.name)")
  println("Local_dir:  $(job.local_dir)")
  println("Server:     $(job.server)")
  println("Server_dir: $(job.server_dir)")
  println("--------------------")
  println("$(length(job.calculations)) calculations:")
  println("")
  for calc in job.calculations
    print_info(calc)
    println("")
  end
end

"""
    print_flags(job::DFJob)

Prints flags of all the calculations in the job.
"""
function print_flags(job::DFJob)
  for calc in job.calculations
    print_flags(calc)
    println("")
  end
end

"""
    print_flags(job::DFJob, calc_filename::String)

Prints flags of the specified calculation.
"""
function print_flags(job::DFJob, calc_filename::String)
  for calc in get_inputs(job,calc_filename)
    print_flags(calc)
    println("")
  end
end

"""
    print_flags(job::DFJob, calc_filenames::Array{String,1})

Prints the flags of the specified calculations.
"""
function print_flags(job::DFJob, calc_filenames::Array{String,1})
  for file in calc_filenames
    print_flags(job,file)
  end
end

"""
    print_flags(job::DFJob, flags)

Prints the specified flags running through all the claculations in the job.
"""
function print_flags(job::DFJob, flags::Array{Symbol,1})
  for flag in flags
    print_flag(job,flag)
  end
end

"""
    print_flag(job::DFJob, flag::Symbol)

Prints the specified flag running through all the calculations in the job.
"""
function print_flag(job::DFJob, flag::Symbol)
  for calc in job.calculations
    print_flag(calc,flag)
  end
end

#all get_inputs return arrays, get_input returns the first element if multiple are found
"""
    get_inputs(job::DFJob, filenames::Array)

Returns an array of the inputs that match one of the filenames.
"""
function get_inputs(job::DFJob, filenames::Array)
  out = DFInput[]
  for name in filenames
    push!(out,filter(x->contains(x.filename,name),job.calculations)...)
  end
  return out
end

"""
    get_inputs(job::DFJob, filename::String)

Returns an array of the input that matches the filename.
"""
function get_inputs(job::DFJob, filename::String)
  return filter(x->contains(x.filename,filename),job.calculations)
end

"""
    get_input(job::DFJob, filename::String)

Returns the input that matches the filename.
"""
function get_input(job::DFJob, filename::String)
  return filter(x->contains(x.filename,filename),job.calculations)[1]
end

"""
    get_input(job::DFJob, filenames::Array)

Returns an array of the inputs that match one of the filenames.
"""
function get_input(job::DFJob, filenames::Array{String,1})
  return get_inputs(job,filenames)
end
#---------------------------------END GENERAL SECTION ------------------#

"""
    change_atoms!(job::DFJob, atoms::Dict{Symbol,<:Array{<:Point3D,1}}, pseudo_set_name=:default, pseudo_fuzzy=nothing)

Sets the data blocks with atomic positions to the new one. This is done for all calculations in the job that have that data. 
If default pseudopotentials are defined, a set can be specified, together with a fuzzy that distinguishes between the possible multiple pseudo strings in the pseudo set. 
These pseudospotentials are then set in all the calculations that need it.
All flags which specify the number of atoms inside the calculation also gets set to the correct value.
"""
function change_atoms!(job::DFJob, atoms::Dict{Symbol,<:Array{<:Point3D,1}};kwargs...)
  for calc in job.calculations
    change_atoms!(calc,atoms;kwargs...)
  end
  nat  = 0
  for pos in values(atoms)
    nat+=length(pos)
  end
  change_flags!(job,:nat=>nat,:ntyp=>length(keys(atoms)))
end

#automatically sets the cell parameters for the entire job, implement others
"""
    change_cell_parameters!(job::DFJob, cell_param::Array{AbstractFloat,2})

Changes the cell parameters in all the input files that have that `DataBlock`.
"""
function change_cell_parameters!(job::DFJob, cell_param::Array{<:AbstractFloat,2})
  for calc in job.calculations
    if typeof(calc) == WannierInput
      alat = get_flag(job,:A)
      change_cell_parameters!(calc,alat*cell_param)
    else
      change_cell_parameters!(calc,cell_param)
    end
  end
end

"""
    change_k_points!(job::DFJob,calc_filename,k_points)

Changes the data in the k point `DataBlock` inside the specified calculation.
"""
function change_k_points!(job::DFJob,calc_filename,k_points)
  change_k_points!(get_input(job,calc_filename),k_points)
end

"""
    change_data_option!(job::DFJob,filenames::Array{String,1}, block_symbol::Symbol,option::Symbol)

Changes the option of specified data block in the specified calculations.
"""
function change_data_option!(job::DFJob,filenames::Array{String,1}, block_symbol::Symbol,option::Symbol)
  for calc in get_inputs(job,filenames)
    change_data_option!(calc,block_symbol,option)
  end
end
change_data_option!(job::DFJob,filename::String, block_symbol::Symbol,option::Symbol) = change_data_option!(job,[filename],block_symbol,option)
"""
    change_data_option!(job::DFJob, block_symbol::Symbol,option::Symbol)

Changes the option of specified data block in all calculations that have the block.
"""
change_data_option!(job::DFJob, block_symbol::Symbol,option::Symbol) = change_data_option!(job,[i.filename for i in job.calculations],block_symbol,option)

"""
    replace_header_word!(job::DFJob,word::String,new_word::String)


Replaces the specified word in the header with the new word.
"""
function replace_header_word!(job::DFJob,word::String,new_word::String)
  for (i,line) in enumerate(job.header)
    if contains(line,word)
      println("Old line:")
      println("   $line")
      job.header[i] = replace(line,word,new_word)
      println("New line:")
      println("   $(job.header[i])")
    end
  end
end
