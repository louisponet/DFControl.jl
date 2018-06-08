"""
    pullfile(server::String, server_dir::String, local_dir::String, filename::String)

Pulls a file from the specified server and server dir to the local dir.
"""
function pullfile(server::String, server_dir::String, local_dir::String, filename::String)
    run(`scp $(server * ":" * server_dir * filename) $local_dir`)
    pulled_files = searchdir(local_dir, filename)
    if !isempty(pulled_files)
        return pulled_files[1]
    end
end

"""
    pullfile(server_dir::String, local_dir::String, filename::String; server=getdefault_server())

Pulls a file from the default server if the default server is specified.
"""
function pullfile(server_dir::String, local_dir::String, filename::String; server=getdefault_server())
    if server != ""
        run(`scp $(server * ":" * server_dir * filename) $local_dir`)
        pulled_files = searchdir(local_dir, filename)
        if !isempty(pulled_files)
            return pulled_files[1]
        end
    else
        error("Define a default server first using 'setdefault_server!(...)'.\n
        Or use function 'pullfile(server::String, server_dir::String, local_dir::String, filename::String)'.")
    end
end

function pullfiles(server_dir::String, local_dir::String, filenames::Vector{String}; serv=getdefault_server())
    pulled_files = String[]
    for file in filenames
       push!(pulled_files, pullfile(server_dir, local_dir, file, serv))
    end
    return pulled_files
end

"""
    pullfile(filepath::String, local_dir::String; server=getdefault_server(), local_filename=nothing)

Pulls a file from the default server if the default server is specified.
"""
function pullfile(filepath::String, local_dir::String; server=getdefault_server(), local_filename=nothing)
    local_dir = form_directory(local_dir)
    if server != ""
        if local_filename != nothing
            run(`scp $(server * ":" * filepath) $(local_dir * local_filename)`)
            return searchdir(local_dir, local_filename)[1]
        else
            run(`scp $(server * ":" * filepath) $(local_dir)`)
            return searchdir(local_dir, splitdir(filepath)[2])[1]
        end
    else
        error("Define a default server first using 'setdefault_server!(...)'.\n
        Or use function 'pullfile(server::String, server_dir::String, local_dir::String, filename::String)'.")
    end
end


#TODO: doesn't work for abinit
#TODO: local pulloutputs!
"""
    pulloutputs(job::DFJob, server="", server_dir="", local_dir=""; job_fuzzy="*job*", extras=String[])

First pulls the job file (specified by job_fuzzy), reads the possible output files and tries to pull them.
Extra files to pull can be specified by the `extras` keyword, works with fuzzy filenames.
"""
function pulloutputs(job::DFJob, server="", server_dir="", local_dir=""; job_fuzzy="*job*", extras=String[], ignore=String[], overwrite=true)
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
    pull_server_file(filename) = pullfile(job.server, job.server_dir, job.local_dir, filename)

    # pull_server_file(job_fuzzy)

    # job_file = searchdir(job.local_dir, strip(job_fuzzy, '*'))[1]
    # inputs, pulloutputs= read_job_filenames(job.local_dir * job_file)
    pulled_outputs = String[]
    for calc in job.inputs
        ofile = outfile(calc)
        if any(contains.(ofile, ignore)) || (ispath(outpath(job,calc)) && !overwrite)
            continue
        end
        try
            pull_server_file(ofile)
            push!(pulled_outputs, ofile)
        end
    end

    for fuzzy in extras
        pull_server_file(fuzzy)
        push!(pulled_outputs, searchdir(job.local_dir, strip(fuzzy,'*'))...)
    end
    return joinpath.(job.local_dir, pulled_outputs)
end

sshcmd(server, cmd) = run(`ssh -t $server $cmd`)
sshreadstring(server, cmd) = readstring(`ssh -t $server $cmd`)
"""
    qstat(server)

If sbatch is running on server it will return the output of the `qstat` command.
"""
qstat(server) = server=="localhost" ? run(`qstat`) : sshcmd(server, "qstat")
qstat()       = qstat(getdefault_server())

"""
    watch_qstat(server)
If sbatch is running on server it will return the output of the `watch qstat` command.
"""
watch_qstat(server) = server=="localhost" ? run(`watch qstat`) : sshcmd(server, "watch qstat")
watch_qstat()       = watch_qstat(getdefault_server())

"Deletes the job from the local or server queue."
qdel(server::String, id::Int) = sshcmd(server, "qdel $id")
qdel(id::Int) = run(`qdel $id`)
qdel(job::DFJob) = runslocal(job) ? qdel(job.metadata[:slurmid]) : qdel(job.server, job.metadata[:slurmid])

function qsub(job::DFJob)
    outstr = ""
    if !runslocal(job)
        push(job)
        outstr = readstring(`ssh -t $(job.server) cd $(job.server_dir) '&&' qsub job.tt`)
    else
        try
            curdir = pwd()
            cd(job.local_dir)
            outstr = readstring(`qsub job.tt`)
            cd(curdir)
        catch
            error("Tried submitting on the local machine but got an error executing `qsub`.")
        end
    end
    return parse(Int, chomp(outstr))
end

"Tests whether a directory exists on a server and if not, creates it."
function mkserverdir(server, dir)
    testfile = joinpath(dir, "tmp.in")
    try
        run(`ssh -t $server touch $testfile`)
        run(`ssh -t $server rm $testfile`)
    catch
        run(`ssh -t $server mkdir -p $dir`)
        info("$dir did not exist on $server, it was created.")
    end
end

"""
    push(job::DFJob)

Pushes a DFJob from it's local directory to its server side directory.
"""
function push(job::DFJob)
    mkserverdir(job.server, job.server_dir)
    scp(file) = run(`scp $(job.local_dir * file) $(job.server * ":" * job.server_dir)`)
    for i in inputs(job)
        scp(i.filename)
    end
    scp("job.tt")
end
