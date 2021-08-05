module Servers
using StructTypes, JLD2, Distributed, HTTP, REPL.TerminalMenus, Parameters
using ..DFControl
using ..Utils

export Server

const SERVER_DIR = DFControl.config_path("servers")

@enum Scheduler Slurm=1 Bash=2

@with_kw mutable struct Server
    name::String           = "nothing"
    username::String       = "nothing"
    domain::String         = "nothing"
    port::Int              = 8080
    scheduler::Scheduler   = Bash
    mountpoint::String     = ""
    julia_exec::String     = "julia"
    default_jobdir::String = homedir()
end

function Server(s::String)
    server = load_server(s) #First check if previous server exists
    if server !== nothing
        return server
    end
    #TODO load existing config 
    # Create new server 
    if occursin("@", s)
        username, domain = split(s, "@")
        name = ""
    else
        @info "Server with name $s not found."
        name = s
        username, domain = "", ""
    end
    @info "Creating new Server configuration..."
    while isempty(name)
        print("Please specify the Server's identifying name:")
        name = readline()
    end
    if load_server(name) !== nothing
        @warn "A server with $name was already configured and will be overwritten."
    end
    while isempty(username)
        print("Username:")
        username = readline()
    end
    while isempty(domain)
        print("Domain:")
        domain = readline()
    end
    print("Port (default: 8080):")
    port_str = readline()
    port = isempty(port_str) ? 8080 : parse(Int, port_str)

    print("Julia Exec (default: julia):")
    julia_str = readline()
    julia = isempty(julia_str) ? "julia" : julia_str

    dir = remotecall_fetch(homedir,
                           Distributed.addprocs([("$username" * "@" * "$domain", 1)])[1])
    print("Jobs top dir (default: $dir):")
    dir = abspath(readline())

    scheduler_choice = request("Please select scheduler:",
                               RadioMenu([string.(instances(Scheduler))...]))
    scheduler_choice == -1 && return
    scheduler = Scheduler(scheduler_choice)
    mounted_choice = request("Has the server been mounted?", RadioMenu(["yes", "no"]))
    mounted_choice == -1 && return
    if mounted_choice == 1
        print("Please specify mounting point:")
        mountpoint = readline()
    else
        mountpoint = ""
    end
    server = Server(name, username, domain, port, scheduler, mountpoint, julia_str, dir)
    println("Server configured as:")
    println(server)
    save(server)
    return server
end

StructTypes.StructType(::Type{Server}) = StructTypes.Mutable()

Base.joinpath(s::Server, p...) = joinpath(s.default_jobdir, p...)

function known_servers(fuzzy = "")
    if ispath(SERVER_DIR)
        servers = [JLD2.load(joinpath(SERVER_DIR, s))["server"]
                   for s in filter(x -> occursin(fuzzy, x), readdir(SERVER_DIR))]
    else
        servers = Server[]
    end
    return servers
end

function load_server(name::String)
    if occursin("@", name)
        return getfirst(x -> x.name == name, known_servers())
    else
        all = known_servers(name)
        if length(all) > 0
            return all[1]
        else
            return nothing
        end
    end
end

function save(s::Server)
    mkpath(SERVER_DIR)
    if ispath(joinpath(SERVER_DIR, s.name * ".jld2"))
        @info "Updating previously existing configuration for server $s."
    else
        JLD2.save(joinpath(SERVER_DIR, s.name * ".jld2"), "server", s)
    end
end

ssh_string(s::Server) = s.username * "@" * s.domain
http_string(s::Server) = "http://$(s.domain):$(s.port)"

function HTTP.request(method::String, s::Server, url, args...; kwargs...)
    return HTTP.request(method, string(http_string(s), url), args...; kwargs...)
end

for f in (:get, :put, :post, :head)
    str = uppercase(string(f))
    @eval function HTTP.$(f)(s::Server, url, args...; kwargs...)
        return HTTP.request("$($str)", s, url, args...; kwargs...)
    end
end

function Distributed.addprocs(server::Server, nprocs::Int = 1, args...; kwargs...)
    if server.name == "localhost"
        proc = Distributed.addprocs(nprocs, args...; kwargs...)
    else
        proc = Distributed.addprocs([(ssh_string(server), nprocs)], args...; kwargs...)
    end
    return proc
end

function remove_server!(name::String)
    return ispath(joinpath(SERVER_DIR, name * ".jld2")) &&
           rm(joinpath(SERVER_DIR, name * ".jld2"))
end

function start(s::Server)
    @info "Starting:\n$s"
    cmd = Cmd(`$(s.julia_exec) --startup-file=no -t auto -e "using DFControl; DFControl.Resource.run($(s.port))"`;
              detach = true)
    proc = addprocs(s)[1]

    p = remotecall(run, proc, cmd; wait = false)
    function isalive()
        try
            HTTP.get(s, "/server_config")
            return true
        catch
            return false
        end
    end

    while !isalive()
        sleep(1)
    end

    @info "Daemon on Server $(s.name) started, listening on port $(s.port)."
    return rmprocs(proc)
end

function maybe_start_server(s::Server)
    try
        HTTP.get(s, "/server_config")
    catch
        start(s)
    end
    return s
end
maybe_start_server(j) = (s = Server(j); return maybe_start_server(s))

kill_server(s::Server) = HTTP.put(s, "/kill_server")

function restart_server(s::Server)
    kill_server(s)
    return start(s)
end

function maybe_create_localhost()
    t = load_server("localhost")
    if t === nothing
        scheduler = Sys.which("sbatch") === nothing ? Bash : Slurm
        if !haskey(ENV, "DFCONTROL_PORT")
            port = 8080
        else
            port = parse(Int, ENV["DFCONTROL_PORT"])
        end
        if !haskey(ENV, "DFCONTROL_JOBDIR")
            dir = homedir()
        else
            dir = ENV["DFCONTROL_PORT"]
        end
        julia_exec = joinpath(Sys.BINDIR, "julia")
        out = Server("localhost", ENV["USER"], "localhost", port, scheduler, "", julia_exec,
                     dir)
        save(out)
        return out
    else
        return t
    end
end

"""
    pull(server::Server, server_file::String, local_file::String)

Pulls `server_file` from the server the `local_file`.
"""
function pull(server::Server, server_file::String, filename::String)
    if server.name == "localhost"
        cp(server_file, filename; force = true)
    else
        run(`scp $(ssh_string(server) * ":" * server_file) $filename`)
    end
end

"""
    push(local_file::String, server::Server, server_file::String)

Pushes the `local_file` to the `server_file` on the server.
"""
function push(filename::String, server::Server, server_file::String)
    if server.name == "localhost"
        cp(filename, server_file; force = true)
    else
        run(`scp $filename $(ssh_string(server) * ":" * server_file)`)
    end
end

end