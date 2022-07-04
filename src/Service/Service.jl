module Service
using Distributed, Pkg, LoggingExtras, Dates, JLD2, LinearAlgebra, JSON3
using ..DFControl: config_path
using ..Utils
using ..Calculations
using ..Jobs
using ..FileIO
using ..Servers
using ..Structures
import ..Database: save, load, exists
# const DAEMON_CONFIG_PATH = config_path("daemon.jld2")
# delete_daemon_config!() = rm(DAEMON_CONFIG_PATH)

PENDING_JOBS_FILE() = config_path("jobs", "pending.txt")
PENDING_WORKFLOWS_FILE() = config_path("workflows", "pending.txt")
QUEUE_FILE() = config_path("jobs", "queue.json")
SLEEP_TIME = 10.0

function daemon_logger()
    p = config_path("logs/daemon")
    mkpath(p)
    FileLogger(joinpath(p, "daemon.log"); append = false)
end

function server_logger()
    p = config_path("logs/runtimes")
    mkpath(p)
    serverid = length(readdir(p)) + 1
    return FileLogger(config_path(joinpath(p, "$serverid.log")); append = false)
end

function restapi_logger()
    p = config_path("logs/daemon")
    mkpath(p)
    FileLogger(joinpath(p, "restapi.log"); append = false)
end
function job_logger(id::Int)
    p = config_path("logs/jobs")
    mkpath(p)
    FileLogger(joinpath(p, "$id.log"))
end

include("schedulers.jl")
include("running.jl")
include("calculation.jl")
include("pseudos.jl")
include("job.jl")
# include("execs.jl")
# include("fileio.jl")
end
