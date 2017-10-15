#TODO this can easily be generalized!!!
#Incomplete this only reads .tt files!!
"""
Loads a Quantum Espresso job from a directory. If no specific input filenames are supplied it will try to find them from a file with name "job".

Input:    job_name::String,
          df_job_dir::String
Optional: T=Float32 -> Type of parsed floats.
Kwargs:   inputs=nothing -> specific input filenames.
          new_home_dir=nothing -> new home directory to store the job later.
          server="" -> server in host@servername format to later send the job to.
          server_dir="" -> server directory where later the job should get pushed to.
"""
function load_qe_job(job_name::String,df_job_dir::String,T=Float32;inputs=nothing,new_homedir=nothing,server="",server_dir="")
  df_job_dir = form_directory(df_job_dir)
  if inputs==nothing
    inputs = read_qe_inputs_from_job_file(df_job_dir*search_dir(df_job_dir,".tt")[1])
  end

  t_calcs = Dict{String,DFInput}()
  flow = Array{Tuple{String,String},1}()
  for (run_command,file) in inputs
    t_calcs[split(file,"_")[end]] = read_qe_input(df_job_dir*file,T)
    push!(flow,(run_command,split(file,"_")[end]))
  end
  if new_homedir!=nothing
    return DFJob(job_name,t_calcs,flow,new_homedir,server,server_dir)
  else
    return DFJob(job_name,t_calcs,flow,df_job_dir,server,server_dir)
  end    
end

"""
Pulls and loads a Quantum Espresso job.

Input: job_name::String, 
       server::String, -> server in host@servername format.
       server_dir::String, -> server side directory.
       local_dir::String -> local directory (will get created if it doesn't exist).
"""
function load_qe_server_job(job_name::String,server::String,server_dir::String,local_dir::String,args...;inputs=nothing)
  pull_job(server,server_dir,local_dir,inputs=inputs)
  return load_qe_job(job_name,local_dir,args...;inputs=inputs,server=server,server_dir=server_dir)
end

#---------------------------------END QUANTUM ESPRESSO SECTION ------------------#

#---------------------------------BEGINNING GENERAL SECTION ---------------------#

#@TODO should we also create a config file for each job with stuff like server etc? and other config things,
#      which if not supplied could contain the default stuff?
"""
Pulls job from server. If no specific inputs are supplied it pulls all .in and .tt files.

Input:  server::String, -> in host@servername format!
        server_dir::String,
        local_dir::String, -> will create the dir if necessary.
Kwargs: inputs=nothing -> specific input filenames.
"""
function pull_job(server::String, server_dir::String, local_dir::String; inputs=nothing)
  server_dir = form_directory(server_dir)
  local_dir  = form_directory(local_dir)
  if !ispath(local_dir)
    mkdir(local_dir)
  end
  if inputs == nothing
    run(`scp -r $(server*":"*server_dir*"*.in") $local_dir`)
    run(`scp -r $(server*":"*server_dir*"*.tt") $local_dir`)
  else
    for input in inputs
      run(`scp -r $(server*":"*server_dir*"*$input") $local_dir`)
    end
  end
end

"""
Saves a DFJob, it's job file and all it's input files.

Input: df_job::DFJob
"""
function save_job(df_job::DFJob)
  if df_job.home_dir == ""
    error("Please specify a valid home_dir!")
  end
  if !ispath(df_job.home_dir)
    mkpath(df_job.home_dir)
  end
  write_job_files(df_job)  
end

#@Incomplete everything is hardcoded for now still!!!! make it configurable
"""
Pushes a DFJob from it's local directory to its server side directory. 

Input: df_job::DFJob
"""
function push_job(df_job::DFJob)
  if !ispath(df_job.home_dir)
    error("Please save the job locally first using save_job(job)!")
  else
    files_to_push = [search_dir(df_job.home_dir,".in");search_dir(df_job.home_dir,".tt")]
    for file in files_to_push
      run(`scp $(df_job.home_dir*file) $(df_job.server*":"*df_job.server_dir)`)
    end
  end
end

#@TODO only uses qsub for now. how to make it more general?
"""
Submit a DFJob. First saves it locally, pushes it to the server then runs the job file on the server.

Input: df_job::DFJob
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
    for calc in values(df_job.calculations)
      for name in fieldnames(calc)[2:end]
        data_dict = getfield(calc,name)
        if typeof(data_dict)<:Dict{Symbol,Dict{Symbol,Any}}
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
Mutatatively change data that is tied to a DFJob. This means that it will run through all the DFInputs and their fieldnames and their Dicts. If it finds a Symbol in one of those that matches a symbol in the new data, it will replace the value of the first symbol with the new value.

Input: df_job::DFJob,
       new_data::Dict{Symbol,Any}
"""
function change_job_data!(df_job::DFJob,new_data::Dict{Symbol,Any})
  for (key,calculation) in df_job.calculations
    for name in fieldnames(calculation)[2:end]
      data_dict = getfield(calculation,name)
      if typeof(data_dict)<:Dict{Symbol,Dict{Symbol,Any}}
        for (block_key,block) in data_dict
          for (flag,value) in block
            if haskey(new_data,flag)
              old_data = value
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
end

#@Incomplete this now assumes that there is only one calculation, would be better to interface with the flow of the DFJob
"""
Sets mutatatively the job data in a calculation block of a DFJob. It will merge the supplied data with the previously present one in the block, changing all previous values to the new ones and adding non-existent ones. 

Input: df_job::DFJob,
       calculation::String, -> calculation in the DFJob.
       block_symbol::Symbol, -> Symbol of the datablock inside the calculation's input file.
       data::Dict{Symbol,Any} -> flags and values to be set.
"""
function set_job_data!(df_job,calculation,block_symbol,data)
  if haskey(df_job.calculations,calculation)
    setfield!(df_job.calculations[calculation],block_symbol,merge(getfield(df_job.calculations[calculation],block_symbol),data))
    println("New input of $block_symbol in $calculation is:\n")
    display(getfield(df_job.calculations[calculation],block_symbol))
    println("\n")
  end
end

"""
Same as above but for multiple calculations.
"""
function set_job_data!(df_job,calculations::Array,block_symbol,data)
  for calculation in calculations
    set_job_data!(df_job,calculation,block_symbol,data)
  end
end
#---------------------------------END GENERAL SECTION ------------------#
