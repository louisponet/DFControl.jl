
#---------------------------------BEGINNING GENERAL SECTION ---------------------#

function pull_file(server::String,server_dir::String,local_dir::String,filename::String)
  run(`scp $(server*":"*server_dir*filename) $local_dir`)
end
"""
    load_job(job_name::String, job_dir::String, T=Float32; job_fuzzy = "job", new_homedir=nothing, server="",server_dir="")

Loads and returns a DFJob. If home_dir is not specified the job directory will ge registered as the local one.
"""
#should we call this load local job?
function load_job(job_name::String, job_dir::String, T=Float32; job_fuzzy = "job", new_homedir=nothing, server="",server_dir="")
  job_dir = form_directory(job_dir)

  t_filenames,t_run_commands,t_should_run = read_job_file(job_dir*search_dir(job_dir,job_fuzzy)[1])
  filenames    = String[]
  run_commands = String[]
  should_run   = Bool[]

  for (i,file) in enumerate(t_filenames)
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
    inputs, run_commands , _ = read_job_file(local_dir*job_file)
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
    load_server_job(job_name::String,server::String,server_dir::String,local_dir::String;job_fuzzy="*job*")

Pulls a server job to local directory and then loads it. A fuzzy search for the job file will be performed and the found input files will be pulled.
"""
function load_server_job(job_name::String, server::String, server_dir::String, local_dir::String; job_fuzzy="*job*")
  pull_job(server,server_dir,local_dir)
  return load_job(job_name,local_dir,server=server,server_dir=server_dir)
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

"""
    check_job_data(df_job,data::Array{Symbol,1})

Check the values of certain flags in a given job if they exist.
"""
function check_job_data(df_job,data_keys)
  out_dict = Dict{Symbol,Any}()
  for s in data_keys
    for (meh,calc) in df_job.calculations
      for name in fieldnames(calc)[2:end]
        data_dict = getfield(calc,name)
        if name == :control_blocks
          for (key,block) in data_dict
            for (flag,value) in block
              if flag == s
                out_dict[s] = value
              end
            end
          end
        else
          for (key,value) in data_dict
            if key == s
              out_dict[s] = value
            end
          end
        end
      end
    end
  end
  return out_dict
end

"""
    change_job_data!(df_job::DFJob,new_data::Dict{Symbol,<:Any})

Mutatatively change data that is tied to a DFJob. This means that it will run through all the DFInputs and their fieldnames and their Dicts.
If it finds a Symbol in one of those that matches a symbol in the new data, it will replace the value of the first symbol with the new value.
"""
function change_job_data!(df_job::DFJob,new_data::Dict{Symbol,<:Any})
  found_keys = Symbol[]
  for (key,calculation) in df_job.calculations
    for name in fieldnames(calculation)[2:end]
      data_dict = getfield(calculation,name)
      if name == :control_blocks
        for (block_key,block) in data_dict
          for (flag,value) in block
            if haskey(new_data,flag)
              old_data = value
              if !(flag in found_keys) push!(found_keys,flag) end
              if typeof(old_data) == typeof(new_data[flag])
                block[flag] = new_data[flag]
                println("$key:\n -> $block_key:\n  -> $flag:\n      $old_data changed to: $(new_data[flag])")
              else
                println("$key:\n -> $block_key:\n  -> $flag:\n    type mismatch old:$old_data ($(typeof(old_data))), new: $(new_data[flag]) ($(typeof(new_data[flag])))\n    Change not applied.")
              end
            end
          end
        end
      else
        for (data_key,data_val) in new_data
          if haskey(data_dict,data_key)
          if !(data_key in found_keys) push!(found_keys,data_key) end
            old_data            = data_dict[data_key]
            if typeof(old_data) == typeof(data_val)
              data_dict[data_key] = data_val
              println("$key:\n -> $name:\n  -> $data_key:\n      $old_data changed to $(data_dict[data_key])")
            else
              println("$key:\n -> $name:\n  -> $data_key:\n    type mismatch old:$old_data ($(typeof(old_data))), new: $data_val ($(typeof(data_val)))\n    Change not applied.")
            end
          end
        end
      end
    end
  end
  for key in found_keys
    pop!(new_data,key)
  end
  if 1 < length(keys(new_data))
    println("flags $(String.(collect(keys(new_data)))) were not found in any input file, please set them first!")
  elseif length(keys(new_data)) == 1
    println("flag '$(String(collect(keys(new_data))[1]))' was not found in any input file, please set it first!")
  end
end

#Incomplete this now assumes that there is only one calculation, would be better to interface with the flow of the DFJob
"""
    set_job_data!(df_job::DFJob,calculation::Int,block_symbol::Symbol,data)

Sets mutatatively the job data in a calculation block of a DFJob. It will merge the supplied data with the previously present one in the block,
changing all previous values to the new ones and adding non-existent ones.
"""
# Input: df_job::DFJob,
#        calculation::String, -> calculation in the DFJob.
#        block_symbol::Symbol, -> Symbol of the datablock inside the calculation's input file.
#        data::Dict{Symbol,Any} -> flags and values to be set.
#Incomplete possibly change calculation to a string rather than an integer but for now it's fine
function set_job_data!(df_job::DFJob, calculation::Int, block_symbol::Symbol, data)
  t_calc = df_job.calculations[calculation][2]
  if block_symbol == :control_blocks
    for (block_key,block_dict) in data
      t_calc.control_blocks[block_key] = merge(t_calc.control_blocks[block_key],data[block_key])
      println("New input of block '$(String(block_key))' in '$(String(block_symbol))' of calculation '$calculation' is now:")
      display(t_calc.control_blocks[block_key])
      println("\n")
    end
  else
    setfield!(t_calc,block_symbol,merge(getfield(t_calc,block_symbol),data))
    println("New input of '$block_symbol' in calculation '$calculation' is:\n")
    display(getfield(t_calc,block_symbol))
    println("\n")
  end
end

"""
    set_job_data!(df_job,calculations::Array,block_symbol,data)

Same as above but for multiple calculations.
"""
function set_job_data!(df_job,calculations::Array,block_symbol,data)
  for calculation in calculations
    set_job_data!(df_job,calculation,block_symbol,data)
  end
end

function remove_job_data!(df_job,calculation,field_symbol,data)
  t_calc = df_job.calculations[calculation][2]
  if field_symbol == :control_blocks
    for d in data
      pop!(t_calc.control_blocks[d])
    end
  else
  end
end

#---------------------------------END GENERAL SECTION ------------------#
