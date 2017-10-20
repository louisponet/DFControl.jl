
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

function pull_outputs(df_job::DFJob, server = "", server_dir = "", local_dir =""; job_fuzzy="*job*",output_fuzzy=nothing)
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

  for (run,output) in zip(should_runs,outputs)
    if run pull_server_file(output) end
  end

  for fuzzy in output_fuzzy
    pull_server_file(fuzzy)
  end
  return outputs
end

function qstat(server)
  run(`ssh -t $server qstat`)
end

function watch_qstat(server)
  run(`ssh -t $server watch qstat`)
end