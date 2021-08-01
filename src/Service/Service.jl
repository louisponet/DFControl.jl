module Service
    using ..DFControl
    using Distributed, Pkg, CodeTracking, LoggingExtras, Dates, JLD2, LinearAlgebra
    using ..DFControl: config_path
    using ..Utils
    # const DAEMON_CONFIG_PATH = config_path("daemon.jld2")
    # delete_daemon_config!() = rm(DAEMON_CONFIG_PATH)

    const RUNNING_JOBS_FILE = config_path("running_jobs.txt")
    const PENDING_JOBS_DIR = config_path("pending_jobs")
    const SERVICE_LOG = config_path("daemon.log")
    const SLEEP_TIME = 10.0

    daemon_logger() = FileLogger(SERVICE_LOG; append = true)

    server_config() = DFC.Server("localhost")

    include("running.jl")
    include("calculation.jl")
    include("job.jl")
    include("versioning.jl")
    include("registry.jl")
    include("execs.jl")
    include("fileio.jl")
    include("slurm.jl")
    
    get_job(p::AbstractString) = DFJob(p)
end
