"""
    pull_file(server::String, server_dir::String, local_dir::String, filename::String)

Pulls a file from the specified server and server dir to the local dir.
"""
function pull_file(server::String, server_dir::String, local_dir::String, filename::String)
  run(`scp $(server*":"*server_dir*filename) $local_dir`)
end

"""
    pull_file(server_dir::String, local_dir::String, filename::String; server=get_default_server())

Pulls a file from the default server if the default server is specified.
"""
function pull_file(server_dir::String, local_dir::String, filename::String; server=get_default_server())
  if server != ""
    run(`scp $(server*":"*server_dir*filename) $local_dir`)
  else
    error("Define a default server first using 'set_default_server!(...)'.\n
    Or use function 'pull_file(server::String, server_dir::String, local_dir::String, filename::String)'.")
  end
end

function read_errors(server::String,server_dir::String; error_fuzzies=["CRASH","*.werr"])
  server_dir = form_directory(server_dir)
  tmp_dir = joinpath(@__DIR__,"tmp")
  if !isdir(tmp_dir)
    mkdir(tmp_dir)
  end
  for fuzzy in error_fuzzies
    run(`scp $(server*":"*server_dir*fuzzy) $local_dir`)
  end

  #for now very dumb!
  crash_readlines = Dict{Symbol,Array{String,1}}()
  for fuzzy in fuzzies
    filenames = search_dir(tmp_dir,strip(fuzzy,'*'))
    if length(filenames)==1
      crash_readlines[filename[1]] = readlines(filenames[1])
    end
  end
  return crash_readlines
end

"""
    pull_outputs(df_job::DFJob, server = "", server_dir = "", local_dir =""; job_fuzzy="*job*",extras=String[])

First pulls the job file (specified by job_fuzzy), reads the possible output files and tries to pull them.
Extra files to pull can be specified by the `extras` keyword, works with fuzzy filenames.
"""
function pull_outputs(df_job::DFJob, server = "", server_dir = "", local_dir =""; job_fuzzy="*job*",extras=String[])
  if df_job.server == "" && server == ""
    error("Error: No job server specified. Please specify it first.")
  elseif server != ""
    df_job.server = server
  end
  if df_job.server_dir == "" && server_dir == ""
    error("Error: No job server_dir specified. Please specify it first.")
  elseif server_dir != ""
    df_job.server_dir = server_dir
  end
  if df_job.local_dir == "" && local_dir == ""
    error("Error: No job local/home directory specified. Please specify it first.")
  elseif server != ""
    df_job.local_dir = local_dir
  end

  pull_server_file(filename) = pull_file(df_job.server,df_job.server_dir,df_job.local_dir,filename)

  pull_server_file(job_fuzzy)

  job_file = search_dir(df_job.local_dir,strip(job_fuzzy,'*'))[1]
  job_name,inputs,outputs,run_command,should_runs = read_job_file(df_job.local_dir*job_file)
  pulled_outputs = String[]
  for (run,output) in zip(should_runs,outputs)
    if run 
      pull_server_file(output)
      push!(pulled_outputs,output)
    end
  end

  for fuzzy in extras
    pull_server_file(fuzzy)
    push!(pulled_outputs,search_dir(df_job.local_dir,strip(fuzzy,'*'))...)
  end
  return pulled_outputs
end

"""
    qstat(server)

If sbatch is running on server it will return the output of the `qstat` command.
"""
function qstat(server)
  run(`ssh -t $server qstat`)
end

qstat() = qstat(default_server())

"""
    watch_qstat(server)
If sbatch is running on server it will return the output of the `watch qstat` command.
"""
function watch_qstat(server)
  run(`ssh -t $server watch qstat`)
end

"Runs `watch_qstat(default_server)`."
watch_qstat() = watch_qstat(default_server())