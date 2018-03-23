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
    pull_file(server_dir::String, local_dir::String, filename::String; server=getdefault_server())

Pulls a file from the default server if the default server is specified.
"""
function pull_file(server_dir::String, local_dir::String, filename::String; server=getdefault_server())
    if server != ""
        run(`scp $(server * ":" * server_dir * filename) $local_dir`)
        pulled_files = search_dir(local_dir, filename)
        if !isempty(pulled_files)
            return pulled_files[1]
        end
    else
        error("Define a default server first using 'setdefault_server!(...)'.\n
        Or use function 'pull_file(server::String, server_dir::String, local_dir::String, filename::String)'.")
    end
end

function pull_files(server_dir::String, local_dir::String, filenames::Vector{String}; server=getdefault_server())
    pulled_files = String[]
    for file in filenames
       push!(pulled_files, pull_file(server_dir, local_dir, file; server))
    end
    return pulled_files
end

"""
    pull_file(filepath::String, local_dir::String; server=getdefault_server(), local_filename=nothing)

Pulls a file from the default server if the default server is specified.
"""
function pull_file(filepath::String, local_dir::String; server=getdefault_server(), local_filename=nothing)
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
        error("Define a default server first using 'setdefault_server!(...)'.\n
        Or use function 'pull_file(server::String, server_dir::String, local_dir::String, filename::String)'.")
    end
end


#TODO: doesn't work for abinit
"""
    pull_outputs(job::DFJob, server="", server_dir="", local_dir=""; job_fuzzy="*job*", extras=String[])

First pulls the job file (specified by job_fuzzy), reads the possible output files and tries to pull them.
Extra files to pull can be specified by the `extras` keyword, works with fuzzy filenames.
"""
function pull_outputs(job::DFJob, server="", server_dir="", local_dir=""; job_fuzzy="*job*", extras=String[])
    if job.server == "" && server == ""
        error("Error: No job server specified. Please specify it first.")
    elseif server != ""
        job.server = server
    end
    if job.server_dir == "" && server_dir == ""
        error("Error: No job server_dir specified. Please specify it first.")
    elseif server_dir != ""
        job.server_dir = server_dir
    end
    if job.local_dir == "" && local_dir == ""
        error("Error: No job local/home directory specified. Please specify it first.")
    elseif server != ""
        job.local_dir = local_dir
    end
    if !ispath(job.local_dir)
        mkpath(job.local_dir)
    end
    pull_server_file(filename) = pull_file(job.server, job.server_dir, job.local_dir, filename)

    # pull_server_file(job_fuzzy)

    # job_file = search_dir(job.local_dir, strip(job_fuzzy, '*'))[1]
    # inputs, outputs= read_job_filenames(job.local_dir * job_file)
    pulled_outputs = String[]
    for calc in job.calculations
        if typeof(calc) == QEInput
            ofile = splitext(calc.filename)[1] * ".out"
        elseif typeof(calc) == WannierInput
            ofile = splitext(calc.filename)[1] * ".wout"
        end
        try
            pull_server_file(ofile)
            push!(pulled_outputs, ofile)
        end
    end

    for fuzzy in extras
        pull_server_file(fuzzy)
        push!(pulled_outputs, search_dir(job.local_dir, strip(fuzzy,'*'))...)
    end
    return job.local_dir .* pulled_outputs
end

"""
    qstat(server)

If sbatch is running on server it will return the output of the `qstat` command.
"""
qstat(server) = run(`ssh -t $server qstat`)
qstat()       = qstat(getdefault_server())

"""
    watch_qstat(server)
If sbatch is running on server it will return the output of the `watch qstat` command.
"""
watch_qstat(server) = run(`ssh -t $server watch qstat`)
watch_qstat()       = watch_qstat(getdefault_server())
