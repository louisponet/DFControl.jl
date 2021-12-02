module Resource
# This module handles all the handling, parsing and transforming HTTP commands,
# and brokers between HTTP requests and Service that fullfils them. 
using HTTP, JSON3, Dates, LoggingExtras, Sockets
using ..DFControl, ..Service
using ..Utils
using ..Servers
using ..Jobs
using ..Calculations
using ..Structures
using ..FileIO

const ROUTER = HTTP.Router()

include("job.jl")
include("fileio.jl")

# GENERAL
kill_server(req) = exit()
HTTP.@register(ROUTER, "PUT", "/kill_server", kill_server)

get_server_config(req) = Service.server_config()
HTTP.@register(ROUTER, "GET", "/server_config", get_server_config)

get_ispath(req) = ispath(job_path(req))
HTTP.@register(ROUTER, "GET", "/get_ispath/*", get_ispath)

get_readdir(req) = readdir(job_path(req))
HTTP.@register(ROUTER, "GET", "/readdir/*", get_readdir)

# PSEUDOS

function pseudos(req)
    fuzzy = String(req.body)
    fuzzy = fuzzy == "\"\"" ? "" : fuzzy
    return Service.pseudos(splitpath(req.target)[end], fuzzy)
end
HTTP.@register(ROUTER, "GET", "/pseudos/*", pseudos)

pseudo_sets(req) = Service.pseudo_sets()
HTTP.@register(ROUTER, "GET", "/pseudos", pseudo_sets)

configure_pseudoset(req) = Service.configure_pseudoset(JSON3.read(req.body,String), job_path(req))
HTTP.@register(ROUTER, "POST", "/configure_pseudoset/*", configure_pseudoset)

rm_pseudos!(req) = Service.rm_pseudos!(JSON3.read(req.body, String))
HTTP.@register(ROUTER, "PUT", "/rm_pseudos", rm_pseudos!)
# EXECS
verify_exec(req) = Service.verify_exec(JSON3.read(req.body, Exec))
HTTP.@register(ROUTER, "GET", "/verify_exec", verify_exec)

known_execs(req) = Service.known_execs(splitpath(req.target)[end])
HTTP.@register(ROUTER, "GET", "/known_execs/*", known_execs)

add_environment(req) = Service.add_environment(JSON3.read(req.body,Jobs.Environment), splitpath(req.target)[end])
HTTP.@register(ROUTER, "POST", "/environment/*", add_environment)
get_environment(req) = Service.get_environment(splitpath(req.target)[end])
HTTP.@register(ROUTER, "GET", "/environment/*", get_environment)

rm_environment!(req) = Service.rm_environment!(splitpath(req.target)[end])
HTTP.@register(ROUTER, "PUT", "/environment/*", rm_environment!)

# RUNNING

function requestHandler(req)
    start = Dates.now()
    @info (timestamp = start, event = "ServiceRequestBegin", tid = Threads.threadid(),
           method = req.method, target = req.target)
    local resp
    try
        obj = HTTP.handle(ROUTER, req)
        if obj === nothing
            resp = HTTP.Response(204)
        else
            resp = HTTP.Response(200, JSON3.write(obj))
        end
    catch e
        s = IOBuffer()
        showerror(s, e, catch_backtrace(); backtrace = true)
        errormsg = String(resize!(s.data, s.size))
        @error errormsg
        resp = HTTP.Response(500, errormsg)
    end
    stop = Dates.now()
    @info (timestamp = stop, event = "ServiceRequestEnd", tid = Threads.threadid(),
           method = req.method, target = req.target, duration = Dates.value(stop - start),
           status = resp.status, bodysize = length(resp.body))
    return resp
end

function run()
    cd(Server("localhost").default_jobdir)
    Service.global_logger(Service.daemon_logger())
    port, server = listenany("0.0.0.0", 8080)
    s = Server("localhost")
    s.port = port
    Servers.save(s)    
    @async HTTP.serve(requestHandler, "0.0.0.0", port, server=server)
    return with_logger(Service.daemon_logger()) do
        Service.main_loop()
    end
end
end
