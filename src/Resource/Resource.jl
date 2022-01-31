module Resource
# This module handles all the handling, parsing and transforming HTTP commands,
# and brokers between HTTP requests and Service that fullfils them. 
using HTTP, JSON3, Dates, LoggingExtras, Sockets, UUIDs
using ..DFControl, ..Service
using ..Utils
using ..Servers
using ..Jobs
using ..Calculations
using ..Structures
using ..FileIO

const ROUTER = HTTP.Router()
const CURRENT_SERVER = Ref{Server}()
const USER_UUID = Ref{UUID}()
include("job.jl")
include("fileio.jl")

# GENERAL

function path(req::HTTP.Request)
    p = req.target
    if !isabspath(p)
        p = joinpath(CURRENT_SERVER[].root_jobdir, p)
    end
    id = findnext(isequal('/'), p, 2)
    if length(p) < id + 1
        return ""
    else
        return p[id+1:end]
    end
end

kill_server(req) = exit()
HTTP.@register(ROUTER, "PUT", "/kill_server", kill_server)

get_server_config(req) = Service.server_config()
HTTP.@register(ROUTER, "GET", "/server_config", get_server_config)

Base.ispath(req::HTTP.Request) = ispath(path(req))
HTTP.@register(ROUTER, "GET", "/ispath/*", ispath)

Base.readdir(req::HTTP.Request) = readdir(path(req))
HTTP.@register(ROUTER, "GET", "/readdir/*", readdir)

Base.mtime(req::HTTP.Request) = mtime(path(req))
HTTP.@register(ROUTER, "GET", "/mtime/*", mtime)

function execute_function(req)
    funcstr = Meta.parse(path(req))
    func = eval(funcstr)
    args = []
    for (t, a) in JSON3.read(req.body, Vector)
        typ = Symbol(t)
        eval(:(arg = JSON3.read($a, $typ)))
        push!(args, arg)
    end
    return func(args...)
end

HTTP.@register(ROUTER, "GET", "/api/*", execute_function)
# PSEUDOS

function pseudos(req)
    fuzzy = String(req.body)
    fuzzy = fuzzy == "\"\"" ? "" : fuzzy
    return Service.pseudos(splitpath(req.target)[end], fuzzy)
end
HTTP.@register(ROUTER, "GET", "/pseudos/*", pseudos)

pseudo_sets(req) = Service.pseudo_sets()
HTTP.@register(ROUTER, "GET", "/pseudos", pseudo_sets)

configure_pseudoset(req) = Service.configure_pseudoset(JSON3.read(req.body,String), path(req))
HTTP.@register(ROUTER, "POST", "/configure_pseudoset/*", configure_pseudoset)

rm_pseudos!(req) = Service.rm_pseudos!(JSON3.read(req.body, String))
HTTP.@register(ROUTER, "PUT", "/rm_pseudos", rm_pseudos!)
# EXECS
verify_exec(req) = Service.verify_exec(JSON3.read(req.body, Exec))
HTTP.@register(ROUTER, "GET", "/verify_exec", verify_exec)

known_execs(req) = (d = JSON3.read(req.body); Service.known_execs(d["exec"],d["dir"]))
HTTP.@register(ROUTER, "GET", "/known_execs/", known_execs)

get_exec(req) = Service.load_exec(splitpath(req.target)[end])
HTTP.@register(ROUTER, "GET", "/exec/*", get_exec)

save_exec(req) = Service.save(JSON3.read(req.body, Exec))
HTTP.@register(ROUTER, "POST", "/exec/", save_exec)

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

function AuthHandler(req)
    if HTTP.hasheader(req, "USER-UUID")
        uuid = HTTP.header(req, "USER-UUID")
        if UUID(uuid) == USER_UUID[]
            return requestHandler(req)
        end
    end
    return HTTP.Response(401, "unauthorized")
end     

function run()
    s = Server("localhost")
    CURRENT_SERVER[] = s
    Service.global_logger(Service.daemon_logger())
    port, server = listenany(ip"0.0.0.0", 8080)
    s.port = port
    Servers.save(s)
    USER_UUID[] = UUID(read(DFC.config_path("user_uuid"), String))
    Threads.@spawn HTTP.serve(AuthHandler, "0.0.0.0", port, server=server)
    with_logger(Service.daemon_logger()) do
        Service.main_loop(s)
    end
    close(server)
    return
    end
end
