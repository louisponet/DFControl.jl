#TODO: doesn't work for abinit
#TODO: local pulloutputs!
"""
    pulloutputs(job::DFJob, server="", server_dir="", dir=""; job_fuzzy="*job*", extras=String[])

First pulls the job file (specified by job_fuzzy), reads the possible output files and tries to pull them.
Extra files to pull can be specified by the `extras` keyword, works with fuzzy filenames.
"""
function pulloutputs(job::DFJob, server = "", server_dir = "", dir = "";
                     job_fuzzy = "*job*", extras = String[], ignore = String[],
                     overwrite = true)
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
    if job.dir == "" && dir == ""
        error("Error: No job local/home directory specified. Please specify it first.")
    elseif server != ""
        job.dir = dir
    end
    if !ispath(job.dir)
        mkpath(job.dir)
    end
    function pull_server_file(filename)
        return pullfile(job.server, job.server_dir, job.dir, filename)
    end

    # pull_server_file(job_fuzzy)

    # job_file = searchdir(job.dir, strip(job_fuzzy, '*'))[1]
    # calculations, pulloutputs= read_job_filenames(job.dir * job_file)
    pulled_outputs = String[]
    for calc in job.calculations
        ofile = outfilename(calc)
        if any(occursin.(ignore, ofile)) || (ispath(outpath(calc)) && !overwrite)
            continue
        end
        try
            pull_server_file(ofile)
            push!(pulled_outputs, ofile)
        catch
            nothing
        end
    end

    for fuzzy in extras
        pull_server_file(fuzzy)
        push!(pulled_outputs, searchdir(job.dir, strip(fuzzy, '*'))...)
    end
    return joinpath.(job.dir, pulled_outputs)
end

sshcmd(server, cmd) = run(`ssh -t $server $cmd`)
sshreadstring(server, cmd) = read(`ssh -t $server $cmd`, String)

"Deletes the job from the local or server queue."
qdel(server::String, id::Int) = sshcmd(server, "qdel $id")
qdel(id::Int) = run(`qdel $id`)
function qdel(job::DFJob)
    if runslocal(job)
        qdel(slurm_jobid(job))
    else
        qdel(job.server, slurm_jobid(job))
    end
end


qsub(job::DFJob) = schedule_job(job, "qsub")
sbatch(job::DFJob) = schedule_job(job, "sbatch")
Base.run(job) = schedule_job(job, "bash")

"Tests whether a directory exists on a server and if not, creates it."
function mkserverdir(server, dir)
    testfile = joinpath(dir, "tmp.in")
    try
        run(`ssh -t $server touch $testfile`)
        run(`ssh -t $server rm $testfile`)
    catch
        run(`ssh -t $server mkdir -p $dir`)
        @info "$dir did not exist on $server, it was created."
    end
end

"""
    pulljob(server::String, server_dir::String, dir::String; job_fuzzy="*job*")

Pulls job from server. If no specific calculations are supplied it pulls all .in and .tt files.
"""
function pulljob(server::String, server_dir::String, dir::String; job_fuzzy = "*job*")
    server_dir = server_dir
    dir  = dir
    if !ispath(dir)
        mkpath(dir)
    end

    pull_server_file(filename) = pullfile(server, server_dir, dir, filename)
    pull_server_file(job_fuzzy)
    job_file = searchdir(dir, strip(job_fuzzy, '*'))[1]

    if job_file != nothing
        calculation_files, output_files = read_job_filenames(joinpath(dir, job_file))
        for file in calculation_files
            pull_server_file(file)
        end
    end
end

pulljob(args...; kwargs...) = pulljob(getdefault_server(), args..., kwargs...)

function push(job::DFJob, s::Server)
    run(`scp $(joinpath(job, "job.jld2")) $(ssh_string(s) * ":" * ".julia/config/DFControl/pending_jobs/" * replace(job.server_dir, "/" => "_")*".jld2")`)
end

Base.joinpath(s::Server, p...) = joinpath(s.default_jobdir, p...)
#Gives the reverse (last job is listed first) of the output, omitting the header lines

