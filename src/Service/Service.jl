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

const RUNNING_JOBS_FILE = config_path("jobs", "running.txt")
const PENDING_JOBS_FILE = config_path("jobs", "pending.txt")
const PENDING_WORKFLOWS_FILE = config_path("workflows", "pending.txt")
const RUNNING_WORKFLOWS_FILE = config_path("workflows", "running.txt")
const SLEEP_TIME = 10.0

daemon_logger() = FileLogger(config_path("logs/daemon.log"); append = false)

function server_logger()
    serverid = length(readdir(config_path("logs/servers"))) + 1
    return FileLogger(config_path("logs/servers/$serverid.log"); append = false)
end

restapi_logger() = FileLogger(config_path("logs/restapi.log"); append = false)
job_logger(id::Int) = FileLogger(config_path("logs/jobs/$id.log"))
next_jobid() = length(readdir(config_path("logs/jobs"))) + 1

server_config() = Server("localhost")
local_server() = server_config()

include("running.jl")
include("calculation.jl")
include("pseudos.jl")
include("job.jl")
# include("execs.jl")
# include("fileio.jl")
include("slurm.jl")
include("bash.jl")
end
