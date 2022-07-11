module Resource
# This module handles all the handling, parsing and transforming HTTP commands,
# and brokers between HTTP requests and Service that fullfils them. 
using HTTP, JSON3, Dates, LoggingExtras, Sockets, UUIDs, ThreadPools
using ..DFControl: config_path
using ..Service
using ..Utils
using ..Servers
using ..Jobs
using ..Calculations
using ..Structures
using ..FileIO
using ..Database

const ROUTER = HTTP.Router()
const CURRENT_SERVER = Ref{Server}()
const USER_UUID = Ref{UUID}()
include("job.jl")
include("fileio.jl")
include("database.jl")

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

get_server_config(req) = Servers.local_server()
HTTP.@register(ROUTER, "GET", "/server_config", get_server_config)

HTTP.@register(ROUTER, "GET", "/isalive", (res) -> true)

Base.ispath(req::HTTP.Request) = ispath(path(req))
HTTP.@register(ROUTER, "GET", "/ispath/*", ispath)

Base.read(req::HTTP.Request) = read(path(req))
HTTP.@register(ROUTER, "GET", "/read/*", read)

Base.write(req::HTTP.Request) = write(path(req), JSON3.read(req.body, Vector{UInt8}))
HTTP.@register(ROUTER, "POST", "/write/*", write)

Base.rm(req::HTTP.Request) = rm(path(req), recursive=true)
HTTP.@register(ROUTER, "POST", "/rm/*", rm)

Base.symlink(req::HTTP.Request) = symlink(JSON3.read(req.body, Vector{String})...)
HTTP.@register(ROUTER, "POST", "/symlink/", symlink)

Base.readdir(req::HTTP.Request) = readdir(path(req))
HTTP.@register(ROUTER, "GET", "/readdir/*", readdir)

Base.mtime(req::HTTP.Request) = mtime(path(req))
HTTP.@register(ROUTER, "GET", "/mtime/*", mtime)

Base.filesize(req::HTTP.Request) = filesize(path(req))
HTTP.@register(ROUTER, "GET", "/filesize/*", filesize)

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

get_exec(req) = Calculations.load(Exec(splitpath(req.target)[end]))
HTTP.@register(ROUTER, "GET", "/execs/*", get_exec)

save_exec(req) = Calculations.save(JSON3.read(req.body, Exec))
HTTP.@register(ROUTER, "POST", "/execs/", save_exec)

add_environment(req) = Jobs.save(JSON3.read(req.body, Jobs.Environment))
HTTP.@register(ROUTER, "POST", "/environment/", add_environment)

get_environment(req) = Jobs.load(Environment(splitpath(req.target)[end]))
HTTP.@register(ROUTER, "GET", "/environment/*", get_environment)

rm_environment!(req) = rm(Environment(splitpath(req.target)[end]))
HTTP.@register(ROUTER, "PUT", "/environment/*", rm_environment!)


# RUNNING
function requestHandler(req)
    start = Dates.now()
    @info (timestamp = string(start), event = "Begin", tid = Threads.threadid(),
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
    @info (timestamp = string(stop), event = "End", tid = Threads.threadid(),
           method = req.method, target = req.target, duration = Dates.value(stop - start),
           status = resp.status, bodysize = length(resp.body))
    return resp
end

function AuthHandler(req)
    if HTTP.hasheader(req, "USER-UUID")
        uuid = HTTP.header(req, "USER-UUID")
        if UUID(uuid) == USER_UUID[]
            t = ThreadPools.spawnbg() do
                requestHandler(req)
            end
            return fetch(t)
        end
    end
    return HTTP.Response(401, "unauthorized")
end     

function run()
    # initialize_config_dir()
    s = Server(gethostname())
    CURRENT_SERVER[] = s
    port, server = listenany(ip"0.0.0.0", 8080)
    s.port = port
    USER_UUID[] = UUID(s.uuid)
    @tspawnat 2 with_logger(Service.server_logger()) do
        Service.main_loop(s)
    end
    Servers.save(s)
    with_logger(Service.restapi_logger()) do
        @info (timestamp = Dates.now(), username = ENV["USER"], host = gethostname(), pid=getpid(), port=port)
        
        HTTP.serve(AuthHandler, "0.0.0.0", port, server=server)
    end
    close(server)
    return
    end
end
