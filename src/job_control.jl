
#---------------------------------BEGINNING GENERAL SECTION ---------------------#

function pull_file(server::String,server_dir::String,local_dir::String,filename::String)
  run(`scp $(server*":"*server_dir*filename) $local_dir`)
end
"""
load_job(job_dir::String, T=Float32; job_fuzzy = "job", new_job_name=nothing, new_homedir=nothing, server="",server_dir="")

Loads and returns a DFJob. If local_dir is not specified the job directory will ge registered as the local one.
"""
#should we call this load local job?
function load_job(job_dir::String, T=Float32; job_fuzzy = "job", new_job_name=nothing, new_local_dir=nothing, server="",server_dir="")
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
  if new_local_dir != nothing
    return DFJob(job_name,t_calcs,new_local_dir,server,server_dir)
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

function add_calculation!(df_job::DFJob,input::DFInput,run_index::Int=length(df_job.calculations)+1;run_command=input.run_command,filename=input.filename)
  input.filename = filename
  input.run_command = run_command
  insert!(df_job.calculations,run_index,input)
end


# """
#     change_job_data!(df_job::DFJob,new_data::Dict{Symbol,<:Any})
#
# Mutatatively change data that is tied to a DFJob. This means that it will run through all the DFInputs and their fieldnames and their Dicts.
# If it finds a Symbol in one of those that matches a symbol in the new data, it will replace the value of the first symbol with the new value.
# """
function change_flags!(df_job::DFJob, new_flag_data::Dict{Symbol,<:Any})
  found_keys = Symbol[]
  for calculation in df_job.calculations
    t_found_keys = change_flags!(calculation,new_flag_data)
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
end

function change_flags!(df_job::DFJob,calc_filenames, new_flag_data::Dict{Symbol,<:Any})
  found_keys = Symbol[]
  for calc in get_inputs(df_job,calc_filenames)
    t_found_keys = change_flags!(calculation,new_flag_data)
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
end

function get_flag(df_job::DFJob,calc_filenames,flag::Symbol)
  for calc in get_inputs(df_job,calc_filenames)
    return get_flag(calc,flag)
  end
end

function get_flag(df_job::DFJob,flag::Symbol)
  for calc in df_job.calculations
    return get_flag(calc,flag)
    # tflag = get_flag(calc,flag)
    # if tflag != nothing
    #   return tflag 
    # end
  end
end

#TODO Change so calculations also have a name.
#TODO change after implementing k_point change so you don't need to specify all this crap
function get_data(df_job::DFJob,calc_filenames,block_symbol::Symbol)
  for calc in get_inputs(df_job,calc_filenames)
    return get_data(calc,block_symbol)
  end
end

function change_data!(df_job::DFJob,calc_filenames, data_block_name::Symbol, new_block_data)
  for calc in get_inputs(df_job,calc_filenames)
    change_data!(calc,data_block_name,new_block_data)
  end
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

#I'm not sure if this is a good idea. Maybe require explicitely also the input filename
#TODO after making defaults this can be done better
function set_flags!(df_job::DFJob,control_block_name::Symbol,data)
  for calc in df_job.calculations
    if :control_blocks in fieldnames(calc)
      set_flags!(calc,control_block_name,data)
    end
  end
end

function set_flags!(df_job::DFJob,data)
  for calc in df_job.calculations
    if :flags in fieldnames(calc)
      set_flags!(calc,data)
    end
  end
end

function remove_flags!(df_job::DFJob,calc_filenames,flags)
  for calc in get_inputs(df_job,calc_filenames)
    remove_flags!(calc,flags)
  end
end

function remove_flags!(df_job::DFJob,flags)
  for calc in df_job.calculations
    remove_flags!(calc,flags)
  end
end

function set_flow!(df_job::DFJob,should_runs::Array{Bool,1})
  assert(length(should_runs)==length(df_job.calculations))
  for (calc,should_run) in zip(df_job.calculations,should_runs)
    calc.run = should_run
  end
  print_flow(df_job)
end

function change_flow!(df_job::DFJob,should_runs::Union{Dict,Array{Tuple{Int,Bool}}})
  for (index,run) in should_runs
    df_job.calculations[index].run = run
  end
  print_flow(df_job)
end

function change_flow!(df_job::DFJob,filenames,should_run)
  for calc in get_inputs(df_job,filenames)
    calc.run = should_run
  end
end

function change_flow!(df_job::DFJob,should_runs::Union{Dict{String,Bool},Array{Tuple{String,Bool}}})
  for (name,should_run) in should_runs
    change_flow!(df_job,name,should_run)
  end
  print_flow(df_job)
end

function change_run_command!(df_job::DFJob,filenames,run_command)
  for calc in get_inputs(df_job,filenames)
    calc.run_command = run_command
    println("Run command of file '$(calc.filename)' is now: '$(calc.run_command)'")
  end
end

function get_run_command(df_job::DFJob,filename)
  for calc in get_inputs(df_job,filename)
    return calc.run_command
  end
end

function print_run_command(df_job::DFJob,filenames)
  for calc in get_inputs(df_job,filenames)
    println("Run command of file '$(calc.filename)' is: '$(calc.run_command)'.")
    println("")
  end
end

function print_flow(df_job::DFJob)
  for (i,calc) in enumerate(df_job.calculations)
    println("$i: $(calc.filename) => runs: $(calc.run)")
  end
end

function print_block(job::DFJob, block_name::Symbol)
  for calc in job.calculations
    if print_block(calc,block_name) println("") end
  end
end

function print_blocks(job::DFJob,calc_filenames)
  for calc in get_inputs(job,calc_filenames)
    print_blocks(calc)
  end
end

function print_blocks(job::DFJob)
  for calc in job.calculations
    print_blocks(calc)
    println("#------------------------------------#")
  end
end

print_data(job::DFJob) = print_blocks(job)
print_data(job::DFJob, calc_filenames) = print_blocks(job, calc_filenames)
print_data(job::DFJob, block_name::Symbol) = print_block(job, block_name)

function print_info(job::DFJob)
  for calc in job.calculations
    print_info(calc)
    println("")
  end
end

function print_flags(job::DFJob)
  for calc in job.calculations
    print_flags(calc)
    println("")
  end
end

function print_flags(job::DFJob,calc_filename::String)
  for calc in get_inputs(job,calc_filename)
    print_flags(calc)
    println("")
  end
end

function print_flags(job::DFJob,calc_filenames::Array{String,1})
  for file in calc_filenames
    print_flags(job,file)
  end
end

function print_flags(job::DFJob,flags)
  for flag in flags
    print_flags(job,flag)
  end
end

function print_flag(job::DFJob,flag::Symbol)
  for calc in job.calculations
    print_flag(calc,flag)
  end
end

#all get_inputs return arrays, get_input returns the first element if multiple are found

function get_inputs(job::DFJob,filenames::Array)
  out = DFInput[]
  for name in filenames
    push!(out,filter(x->contains(x.filename,name),job.calculations)...)
  end
  return out
end

function get_inputs(job::DFJob,filename::String)
  return filter(x->contains(x.filename,filename),job.calculations)
end

function get_input(job::DFJob,filename::String)
  return filter(x->contains(x.filename,filename),job.calculations)[1]
end

function get_input(job::DFJob,filenames::Array{String,1})
  return get_inputs(job,filenames)
end
#---------------------------------END GENERAL SECTION ------------------#

#TODO after adding defaults can automatically get correct pseudos as well
function change_atoms!(job::DFJob,atoms::Dict{Symbol,<:Array{<:Point3D,1}})
  for calc in job.calculations
    change_atoms!(calc,atoms)
  end
end

#automatically sets the cell parameters for the entire job, implement others
function change_cell_parameters!(job::DFJob,cell_param::Array{AbstractFloat,2})
  for calc in job.calculations
    if typeof(calc) == WannierInput
      alat = get_flag(job,:A)
      change_cell_parameters!(calc,alat*cell_param)
    else
      change_cell_parameters!(calc,cell_param)
    end
  end
end

function change_k_points!(job::DFJob,calc_filename,k_points)
  change_k_points!(get_input(job,calc_filename),k_points)
end