"""
    pull_file(server::String, server_dir::String, local_dir::String, filename::String)

Pulls a file from the specified server and server dir to the local dir.
"""
function pull_file(server::String, server_dir::String, local_dir::String, filename::String)
    run(`scp $(server * ":" * server_dir * filename) $local_dir`)
    pulled_files = search_dir(local_dir, filename)
    if !isempty(pulled_files)
        return pulled_files[1]
    end
end

"""
    pull_file(server_dir::String, local_dir::String, filename::String; server=get_default_server())

Pulls a file from the default server if the default server is specified.
"""
function pull_file(server_dir::String, local_dir::String, filename::String; server=get_default_server())
    if server != ""
        run(`scp $(server * ":" * server_dir * filename) $local_dir`)
        pulled_files = search_dir(local_dir, filename)
        if !isempty(pulled_files)
            return pulled_files[1]
        end
    else
        error("Define a default server first using 'set_default_server!(...)'.\n
        Or use function 'pull_file(server::String, server_dir::String, local_dir::String, filename::String)'.")
    end
end

function pull_files(server_dir::String, local_dir::String, filenames::Array{String}; server=get_default_server())
    pulled_files = String[]
    for file in filenames
       push!(pulled_files, pull_file(server_dir, local_dir, file; server))
    end 
    return pulled_files
end

"""
    pull_file(filepath::String, local_dir::String; server=get_default_server(), local_filename=nothing)

Pulls a file from the default server if the default server is specified.
"""
function pull_file(filepath::String, local_dir::String; server=get_default_server(), local_filename=nothing)
    local_dir = form_directory(local_dir)
    if server != ""
        if local_filename != nothing
            run(`scp $(server * ":" * filepath) $(local_dir * local_filename)`)
            return search_dir(local_dir, local_filename)[1]
        else
            run(`scp $(server * ":" * filepath) $(local_dir)`)
            return search_dir(local_dir, splitdir(filepath)[2])[1]
        end
    else
        error("Define a default server first using 'set_default_server!(...)'.\n
        Or use function 'pull_file(server::String, server_dir::String, local_dir::String, filename::String)'.")
    end
end

"""
    pull_outputs(df_job::DFJob, server="", server_dir="", local_dir=""; job_fuzzy="*job*", extras=String[])

First pulls the job file (specified by job_fuzzy), reads the possible output files and tries to pull them.
Extra files to pull can be specified by the `extras` keyword, works with fuzzy filenames.
"""
function pull_outputs(df_job::DFJob, server="", server_dir="", local_dir=""; job_fuzzy="*job*", extras=String[])
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
    if !ispath(df_job.local_dir)
        mkpath(df_job.local_dir)
    end
    pull_server_file(filename) = pull_file(df_job.server, df_job.server_dir, df_job.local_dir, filename)
    
    pull_server_file(job_fuzzy)
    
    job_file = search_dir(df_job.local_dir, strip(job_fuzzy, '*'))[1]
    job_data = read_job_file(df_job.local_dir * job_file)
    pulled_outputs = String[]
    for (run, output) in zip(job_data[:should_run], job_data[:output_files])
        if run 
            pull_server_file(output)
            push!(pulled_outputs, output)
        end
    end
    
    for fuzzy in extras
        pull_server_file(fuzzy)
        push!(pulled_outputs, search_dir(df_job.local_dir, strip(fuzzy,'*'))...)
    end
    return df_job.local_dir .* pulled_outputs
end

"""
    qstat(server)

If sbatch is running on server it will return the output of the `qstat` command.
"""
function qstat(server)
    run(`ssh -t $server qstat`)
end
qstat() = qstat(get_default_server())

"""
    watch_qstat(server)
If sbatch is running on server it will return the output of the `watch qstat` command.
"""
function watch_qstat(server)
    run(`ssh -t $server watch qstat`)
end
watch_qstat() = watch_qstat(get_default_server())