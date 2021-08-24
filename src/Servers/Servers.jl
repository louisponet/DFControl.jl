module Servers
using StructTypes, Distributed, JSON3, HTTP, REPL.TerminalMenus, Parameters
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
    try
        @info "Trying to pull existing configuratiaon from $username@$domain..."
        localpath = joinpath(SERVER_DIR, name*".json")
        remotepath = ".julia/config/DFControl/servers/localhost.json"
        run(`scp $(username * "@" * domain):$remotepath $localpath`)
        tserver = load_server(name)
        tserver.name = name
        tserver.username= username
        tserver.domain=domain
        println("Server configured as:")
        println(tserver)
        save(tserver)
        return tserver
    catch
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

        @info "Trying to connect to $username@$domain..."
        dir = remotecall_fetch(homedir,
                               Distributed.addprocs([("$username" * "@" * "$domain", 1)], exename=julia)[1])
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
end

StructTypes.StructType(::Type{Server}) = StructTypes.Mutable()

Base.joinpath(s::Server, p...) = joinpath(s.default_jobdir, p...)

function known_servers(fuzzy = "")
    if ispath(SERVER_DIR)
        servers = [JSON3.read(read(joinpath(SERVER_DIR, s), String), Server)
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
    if ispath(joinpath(SERVER_DIR, s.name * ".json"))
        @info "Updating previously existing configuration for server $s."
    end
    JSON3.write(joinpath(SERVER_DIR, s.name * ".json"),  s)
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
    if server.domain == "localhost"
        proc = Distributed.addprocs(nprocs, args...; kwargs...)
    else
        proc = Distributed.addprocs([(ssh_string(server), nprocs)], args...; exename=`$(server.julia_exec)`, dir=server.default_jobdir, tunnel=true, kwargs...)
    end
    return proc
end

function remove_server!(name::String)
    return ispath(joinpath(SERVER_DIR, name * ".json")) &&
           rm(joinpath(SERVER_DIR, name * ".json"))
end

function isalive(s::Server)
    try
        HTTP.get(s, "/server_config")
        return true
    catch
        return false
    end
end

function start(s::Server)
    @info "Starting:\n$s"
    cmd = Cmd(`$(s.julia_exec) --startup-file=no -t auto -e "using DFControl; DFControl.Resource.run($(s.port))"`;
              detach = true)
    proc = addprocs(s)[1]

    p = remotecall(run, proc, cmd; wait = false)

    while !isalive(s)
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
    if server.domain == "localhost"
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
    if server.domain == "localhost"
        cp(filename, server_file; force = true)
    else
        run(`scp $filename $(ssh_string(server) * ":" * server_file)`)
    end
end

function server_command(s::Server, cmd)
    out = Pipe()
    err = Pipe()
    if s.domain == "localhost"
        process = run(pipeline(ignorestatus(cmd), stdout=out, stderr=err))
    else
        process = run(pipeline(ignorestatus(`ssh $(ssh_string(s)) source /etc/profile '&''&' $cmd`), stdout=out, stderr=err))
    end
    close(out.in)
    close(err.in)

    stdout = String(read(out))
    stderr = String(read(err))
    return (
      stdout = stdout,
      stderr = stderr,
      code = process.exitcode
    )            
end

function has_modules(s::Server)
    try 
        server_command(s, `module avail`).code == 0
    catch
        false
    end
end

function available_modules(s::Server)
    if has_modules(s)
        return server_command(s, `module avail`)
    else
        return String[]
    end
end

end
