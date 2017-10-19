
#---------------------------------BEGINNING GENERAL SECTION ---------------------#

function pull_file(server::String,server_dir::String,local_dir::String,filename::String)
  run(`scp $(server*":"*server_dir*filename) $local_dir`)
end
"""
    load_job(job_dir::String, T=Float32; job_fuzzy = "job", new_job_name=nothing, new_homedir=nothing, server="",server_dir="")

Loads and returns a DFJob. If home_dir is not specified the job directory will ge registered as the local one.
"""
#should we call this load local job?
function load_job(job_dir::String, T=Float32; job_fuzzy = "job", new_job_name=nothing, new_homedir=nothing, server="",server_dir="")
  job_dir = form_directory(job_dir)

  job_name,t_inputs,t_outputs,t_run_commands,t_should_run = read_job_file(job_dir*search_dir(job_dir,job_fuzzy)[1])
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
      push!(t_calcs,read_wannier_input(filename*".win",T,run_command=run_command,run=run))
    else
      push!(t_calcs,read_qe_input(filename,T,run_command=run_command,run=run))
    end
  end
  if new_homedir != nothing
    return DFJob(job_name,t_calcs,new_homedir,server,server_dir)
  else
    return DFJob(job_name,t_calcs,job_dir,server,server_dir)
  end
end

#TODO should we also create a config file for each job with stuff like server etc? and other config things,
#      which if not supplied could contain the default stuff?
"""
    pull_job(server::String, server_dir::String, local_dir::String; job_fuzzy="*job*")
Pulls job from server. If no specific inputs are supplied it pulls all .in and .tt files.
"""
# Input:  server::String, -> in host@servername format!
#         server_dir::String,
#         local_dir::String, -> will create the dir if necessary.
# Kwargs: inputs=nothing -> specific input filenames.
function pull_job(server::String, server_dir::String, local_dir::String; job_fuzzy="*job*")
  server_dir = form_directory(server_dir)
  local_dir  = form_directory(local_dir)
  if !ispath(local_dir)
    mkdir(local_dir)
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

"""
    load_server_job(server::String, server_dir::String, local_dir::String; job_fuzzy="*job*", job_name=nothing)

Pulls a server job to local directory and then loads it. A fuzzy search for the job file will be performed and the found input files will be pulled.
"""
function load_server_job(server::String, server_dir::String, local_dir::String; job_fuzzy="*job*", new_job_name=nothing)
  pull_job(server,server_dir,local_dir)
  return load_job(local_dir,server=server,server_dir=server_dir,new_job_name = new_job_name)
end

"""
    save_job(df_job::DFJob)

Saves a DFJob, it's job file and all it's input files.
"""
function save_job(df_job::DFJob)
  home_dir = df_job.home_dir
  if home_dir == ""
    error("Please specify a valid home_dir!")
  end
  home_dir = form_directory(df_job.home_dir)
  if !ispath(home_dir)
    mkpath(home_dir)
  end
  df_job.home_dir = home_dir
  write_job_files(df_job)
end

#Incomplete everything is hardcoded for now still!!!! make it configurable
"""
    push_job(df_job::DFJob)

Pushes a DFJob from it's local directory to its server side directory.
"""
function push_job(df_job::DFJob)
  if !ispath(df_job.home_dir)
    error("Please save the job locally first using save_job(job)!")
  else
    for calc in df_job.calculations
      run(`scp $(df_job.home_dir*calc.filename) $(df_job.server*":"*df_job.server_dir)`)
    end
    run(`scp $(df_job.home_dir*"job.tt") $(df_job.server*":"*df_job.server_dir)`)
  end
end

#TODO only uses qsub for now. how to make it more general?
"""
    submit_job(df_job::DFJob)

Submit a DFJob. First saves it locally, pushes it to the server then runs the job file on the server.
"""
function submit_job(df_job::DFJob)
  if df_job.server == ""
    error("Please specify a valid server address first!")
  elseif df_job.server_dir == ""
    error("Please specify a valid server directory first!")
  end
  save_job(df_job)
  push_job(df_job)
  run(`ssh -t $(df_job.server) cd $(df_job.server_dir) '&&' qsub job.tt`)
end

function add_calculation!(df_job::DFJob,input::DFInput,run_index::Int=length(df_job.calculations)+1)
  insert!(df_job.calculations,run_index,input)
end

# All the methods to change the inpût control flags, if you want to implement another kind of calculation add a similar one here!

function change_input_control_flags!(input::QEInput, new_flag_data::Dict{Symbol,<:Any})
  found_keys = Symbol[]
  for block in input.control_blocks
    for (flag,value) in new_flag_data
      if haskey(block.flags,flag)
        old_data = block.flags[flag]
        if !(flag in found_keys) push!(found_keys,flag) end
        if typeof(block.flags[flag]) == typeof(new_flag_data[flag])
          block.flags[flag] = new_flag_data[flag]
          println("$(input.filename):\n -> $(block.name):\n  -> $flag:\n      $old_data changed to: $(new_flag_data[flag])")
        else
          println("$(input.filename):\n -> $(block.name):\n  -> $flag:\n     type mismatch old:$old_data ($(typeof(old_data))), new: $(new_flag_data[flag]) ($(typeof(new_flag_data[flag])))\n    Change not applied.")
        end
      end
    end
  end
  return found_keys
end

function change_input_control_flags!(input::WannierInput, new_flag_data::Dict{Symbol,<:Any})
  found_keys = Symbol[]
  for (flag,value) in new_flag_data
    if haskey(input.flags,flag)
      old_data = input.flags[flag]
      if !(flag in found_keys) push!(found_keys,flag) end
      if typeof(input.flags[flag]) == typeof(new_flag_data[flag])
        input.flags[flag] = new_flag_data[flag]
        println("$(input.filename):\n -> $flag:\n      $old_data changed to: $(new_flag_data[flag])")
      else
        println("$(input.filename):\n -> $flag:\n     type mismatch old:$old_data ($(typeof(old_data))), new: $(new_flag_data[flag]) ($(typeof(new_flag_data[flag])))\n    Change not applied.")
      end
    end
  end
  return found_keys
end

function change_input_data!(input::DFInput, block_name::Symbol, new_block_data)
  for data_block in input.data_blocks
    if data_block.name == block_name
      if typeof(data_block.data) == typeof(new_block_data)
        old_data = data_block.data
        data_block.data = new_block_data
        println("Block data '$(data_block.name)' in input  '$(input.filename)' is now:\n")
        display(data_block.data)
      end
    end
  end
end
# """
#     change_job_data!(df_job::DFJob,new_data::Dict{Symbol,<:Any})
#
# Mutatatively change data that is tied to a DFJob. This means that it will run through all the DFInputs and their fieldnames and their Dicts.
# If it finds a Symbol in one of those that matches a symbol in the new data, it will replace the value of the first symbol with the new value.
# """
function change_job_control_flags!(df_job::DFJob, new_flag_data::Dict{Symbol,<:Any})
  found_keys = Symbol[]
  for calculation in df_job.calculations
    t_found_keys = change_input_control_flags!(calculation,new_flag_data)
    for key in t_found_keys
      if !(key in found_keys) push!(found_keys,key) end
    end
  end
  for key in found_keys
    pop!(new_flag_data,key)
  end
  if 1 < length(keys(new_flag_data))
    println("flags $(String.(collect(keys(new_flag_data)))) were not found in any input file, please set them first!")
  elseif length(keys(new_flag_data)) == 1
    println("flag '$(String(collect(keys(new_flag_data))[1]))' was not found in any input file, please set it first!")
  end
  return df_job
end

function change_job_data!(df_job::DFJob,calculations::Array{Int,1} data_block_name::Symbol, new_block_data)
  for calc in calculations
    change_input_data!(df_job.calculations[calc],data_block_name,new_block_data)
  end
  return df_job
end

#here comes the code for all the setting of flags of different inputs
function set_input_control_flags!(input::QEInput,control_block_name::Symbol,data)
  for block in input.control_blocks
    if block.name == control_block_name
      block.flags = merge((x,y) -> typeof(x) == typeof(y) ? y : x,block.flags,data)
      println("New input of block '$(block.name)' of calculation '$(input.filename)' is now:")
      display(block.flags)
      println("\n")
    end
  end
end

function set_input_control_flags!(input::WannierInput,data)
  input.flags = merge((x,y) -> typeof(x) == typeof(y) ? y : x,block.flags,data)
  println("New input of calculation '$(input.filename)' is now:")
  display(input.flags)
  println("\n")
end
#Incomplete this now assumes that there is only one calculation, would be better to interface with the flow of the DFJob
# """
#     set_job_data!(df_job::DFJob,calculation::Int,block_symbol::Symbol,data)
#
# Sets mutatatively the job data in a calculation block of a DFJob. It will merge the supplied data with the previously present one in the block,
# changing all previous values to the new ones and adding non-existent ones.
# """
# Input: df_job::DFJob,
#        calculation::String, -> calculation in the DFJob.
#        block_symbol::Symbol, -> Symbol of the datablock inside the calculation's input file.
#        data::Dict{Symbol,Any} -> flags and values to be set.
#Incomplete possibly change calculation to a string rather than an integer but for now it's fine
function set_job_control_flags!(df_job::DFJob,control_block_name::Symbol,data)
  for calc in df_job.calculations
    if :control_blocks in fieldnames(calc) && control_block_name != nothing
      set_input_control_flags!(calc,control_block_name,data)
    end
  end
  return df_job
end

function set_job_control_flags!(df_job::DFJob,data)
  for calc in df_job.calculations
    if :flags in fieldnames(calc)
      set_input_control_flags!(calc,data)
    end
  end
  return df_job
end

#removes an input control flag, if you want to implement another input add a similar function here!
function remove_input_control_flag!(input::QEInput,flag)
  for block in input.control_blocks
    if haskey(block.flags,flag)
      pop!(block.flags,flag)
      println("Removed flag '$flag' from block '$(block.name)' in input '$(input.filename)'")
    end
  end
end

function remove_input_control_flag!(input::WannierInput,flag)
  if haskey(input.flags,flag)
    pop!(input.flags,flag,false)
    println("Removed flag '$flag' from input '$(input.filename)'")
  end
end

function remove_job_control_flag!(df_job::DFJob,calculation,flag)
  remove_input_control_flag!(df_job.calculations[calculation],flag)
  return df_job
end

function remove_job_control_flag!(df_job::DFJob,calculations::Array,flag)
  for calc in calculations
    remove_job_control_flag!(df_job,calc,flag)
  end
  return df_job
end

function remove_job_control_flags!(df_job::DFJob, calculation::Int, flags::Array{Symbol,1})
  remove_job_control_flag!.(df_job,calculation,flags)
  return df_job
end

function remove_job_control_flags!(df_job::DFJob, calculations::Array{Int,1}, flags::Array{Symbol,1})
  for calculation in calculations
    remove_job_control_flag!.(df_job,calculation,flags)
  end
  return df_job
end

function remove_job_control_flag!(df_job::DFJob,flag)
  for calculation in df_job.calculation
    remove_input_control_flag!(calc,flag)
  end
  return df_job
end

function remove_job_control_flags!(df_job::DFJob,flags::Array{Symbol,1})
  for flag in flags
    for calculation in df_job.calculations
      remove_input_control_flag!(calculation,flag)
    end
  end
  return df_job
end


function set_should_run!(df_job::DFJob,should_runs::Array{Bool,1})
  assert(length(should_runs)==length(df_job.calculations))
  for (calc,should_run) in zip(df_job.calculations,should_runs)
    calc.run = should_run
  end
  return df_job
end
#---------------------------------END GENERAL SECTION ------------------#
