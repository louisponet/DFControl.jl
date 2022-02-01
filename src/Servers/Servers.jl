module Servers
using StructTypes, Distributed, JSON3, HTTP, REPL.TerminalMenus, Parameters, UUIDs
using ..DFControl: config_path
using ..Utils
using ..Jobs

export Server

const SERVER_DIR = config_path("servers")

@enum Scheduler Slurm=1 Bash=2

"""
    Server(name::String, username::String, domain::String, port::Int, scheduler::Scheduler, mountpoint::String,
           julia_exec::String, root_jobdir::String, local_port::Int, max_concurrent_jobs::Int)
    Server(name::String)

A [`Server`](@ref) represents a remote daemon that has the label `name`. It runs on the server defined by
`username` and `domain`. The requirement is that `ssh` is set up in such a way that `ssh username@domain` is
possible, i.e. ssh-copy-id must have been used to not require passwords while executing `ssh` commands.

The daemon will listen to the `port` for http requests and if `local_port` is specified,
a tunnel will be created to guarantee a connection. This is useful in the case that the login node on the remote
server can change.

Calling [`Server`](@ref) with a single `String` will either load the configuration that was previously saved with that label, or go through an interactive setup of a new server.
"""
mutable struct Server
    name::String
    username::String
    domain::String
    port::Int
    scheduler::Scheduler
    mountpoint::String
    julia_exec::String
    root_jobdir::String
    local_port::Int
    max_concurrent_jobs::Int
end

function Server(s::String)
    server = load(s) #First check if previous server exists
    if server !== nothing
        return server
    end
    # Create new server 
    @info "Creating new Server configuration..."
    if occursin("@", s)
        username, domain = split(s, "@")
        name = ask_input(String, "Please specify the Server's identifying name:")
    else
        username = ask_input(String, "Username")
        domain = ask_input(String, "Domain")
        name = s
    end
    if load(name) !== nothing
        @warn "A server with $name was already configured and will be overwritten."
    end
    @info "Trying to pull existing configuration from $username@$domain..."

    server = load_remote_config(username, domain)
    if server !== nothing

        local_choice = request("Should a local tunnel be created?", RadioMenu(["yes", "no"]))
        local_choice == -1 && return
        if local_choice == 1
            loc_port = ask_input(Int, "Local port", 8123)
        else
            loc_port = 0
        end
             
        server.name       = name
        server.domain     = domain
        server.local_port = loc_port
    else
        @info "Couldn't pull server configuration, creating new..."
        port  = ask_input(Int, "Port", 8080)
        julia = ask_input(String, "Julia Exec", "julia")
        @info "Trying to connect to $username@$domain..."
        if julia != "julia"
            while server_command(username, domain, `ls $julia`).exitcode != 0
                @warn "$julia, no such file or directory."
                julia = ask_input(String, "Julia Exec")
            end
        end

        hdir = server_command(username, domain, `pwd`).stdout
        dir = ask_input(String, "Default Jobs directory", hdir)
        if dir != hdir
            while server_command(username, domain, `ls $dir`).exitcode != 0
                @warn "$dir, no such file or directory."
                dir = ask_input(String, "Default Jobs directory")
            end
        end

        scheduler_choice = request("Please select scheduler:",
                                   RadioMenu([string.(instances(Scheduler))...]))
                                   
        scheduler_choice == -1 && return
        scheduler = Scheduler(scheduler_choice)
        mounted_choice = request("Has the server been mounted?", RadioMenu(["yes", "no"]))
        mounted_choice == -1 && return
        if mounted_choice == 1
            mountpoint = ask_input(String, "Mounting point")
        else
            mountpoint = ""
        end
        local_choice = request("Should a local tunnel be created?", RadioMenu(["yes", "no"]))
        local_choice == -1 && return
        if local_choice == 1
            loc_port = ask_input(Int, "Local port", 8123)
        else
            loc_port = 0
        end
        max_concurrent_jobs = ask_input(Int, "Max Concurrent Jobs", 100)
        
        server = Server(name, username, domain, port, scheduler, mountpoint, julia, dir, loc_port, max_concurrent_jobs)

    end
    save(server)
    return server
end

StructTypes.StructType(::Type{Server}) = StructTypes.Struct()
islocal(s::Server) = s.domain == "localhost"

Base.joinpath(s::Server, p...) = joinpath(s.root_jobdir, p...)
Base.ispath(s::Server, p...) =
    JSON3.read(HTTP.get(s, "/ispath/" * joinpath(p...)).body, Bool)

Utils.searchdir(s::Server, dir, str) = joinpath.(dir, filter(x->occursin(str, x), readdir(s, dir))) 

parse_config(config) = JSON3.read(config, Server)
read_config(config_file) = parse_config(read(config_file, String))

function load_remote_config(username, domain; name="localhost")
    cmd = "cat ~/.julia/config/DFControl/servers/$name.json"
    t = server_command(username, domain, cmd)
    if t.exitcode != 0
        return nothing
    else
        return parse_config(t.stdout)
    end
end

function known_servers(fuzzy = "")
    if ispath(SERVER_DIR)
        servers = [JSON3.read(read(joinpath(SERVER_DIR, s), String), Server)
                   for s in filter(x -> occursin(fuzzy, x), readdir(SERVER_DIR))]
    else
        servers = Server[]
    end
    return servers
end

function load(name::String)
    if occursin("@", name)
        n = split(name, "@")[end]
        return getfirst(x -> x.name == n, known_servers())
    else
        all = known_servers(name)
        if length(all) > 0
            return all[1]
        else
            return nothing
        end
    end
end

"""
    save(s::Server)

Saves the server configuration for later use. Can be loaded using [`Server`](@ref) with the name of `s`.
If the server is not local it will also upload the configuration to the remote host.
"""
function save(s::Server)
    mkpath(SERVER_DIR)
    if ispath(joinpath(SERVER_DIR, s.name * ".json"))
        @info "Updating previously existing configuration for server $s."
    end
    JSON3.write(joinpath(SERVER_DIR, s.name * ".json"),  s)
    if !islocal(s)
        @info "Uploading server configuration."
        t = deepcopy(s)
        t.domain = "localhost"
        t.local_port = 0
        t.name = "localhost"
        tf = tempname()
        JSON3.write(tf,  t)
        push(tf, s, "~/.julia/config/DFControl/servers/localhost.json")
    end
end

ssh_string(s::Server) = s.username * "@" * s.domain
http_string(s::Server) = s.local_port != 0 ? "http://localhost:$(s.local_port)" : "http://$(s.domain):$(s.port)"

function HTTP.request(method::String, s::Server, url, args...; kwargs...)
    uuid = read(config_path("user_uuid"), String)
    return HTTP.request(method, string(http_string(s), url), ["USER-UUID" => uuid], args...; kwargs...)
end

for f in (:get, :put, :post, :head)
    str = uppercase(string(f))
    @eval function HTTP.$(f)(s::Server, url, args...; kwargs...)
        return HTTP.request("$($str)", s, url, args...; kwargs...)
    end
end

function Distributed.addprocs(server::Server, nprocs::Int = 1, args...; kwargs...)
    if islocal(server)
        proc = Distributed.addprocs(nprocs, args...; kwargs...)
    else
        proc = Distributed.addprocs([(ssh_string(server), nprocs)], args...; exename=`$(server.julia_exec)`, dir=server.root_jobdir, tunnel=true, kwargs...)
    end
    return proc
end

function rm(s::Server)
    return ispath(joinpath(SERVER_DIR, s.name * ".json")) &&
           rm(joinpath(SERVER_DIR, s.name * ".json"))
end

function destroy_tunnel(s)
    t = getfirst(x->occursin("ssh -N -f -L $(s.local_port)", x), split(read(pipeline(`ps aux` , stdout = `grep $(s.local_port)`), String), "\n"))
        
    if t !== nothing
        run(`kill $(split(t)[2])`)
    end
end
function construct_tunnel(s)
    run(Cmd(`ssh -N -f -L $(s.local_port):localhost:$(s.port) $(ssh_string(s))`, detach=true))
end

"""
    isalive(s::Server)

Will try to fetch some data from `s`. If the server is not running this will fail and
the return is `false`.
"""
function isalive(s::Server)
    try
        resp = HTTP.get(s, "/server_config", readtimeout=15, retry=false)
        return JSON3.read(resp.body, Server).username == s.username
    catch
        tserver = load_remote_config(s.username, s.domain)
        s.port = tserver.port
        if s.local_port != 0
            destroy_tunnel(s)
            construct_tunnel(s)
        end
            
        try 
            resp = HTTP.get(s, "/server_config", readtimeout=15, retry=false)
            return JSON3.read(resp.body, Server).username == s.username
        catch
            return false
        end
    end
end

"""
    start(s::Server)

Launches the daemon process on  the host [`Server`](@ref) `s`.
"""
function start(s::Server)
    @assert !isalive(s) "Server is already up and running."
    @info "Starting:\n$s"
    if islocal(s)
        t = ispath(config_path("self_destruct"))
    else
        cmd = "cat ~/.julia/config/DFControl/self_destruct"
        t = server_command(s, cmd).exitcode == 0
    end
        
    @assert !t "Self destruction was previously triggered, signalling issues on the Server.\nPlease investigate and if safe, remove ~/.julia/config/DFControl/self_destruct"

    # Here we clean up previous connections and commands
    if !islocal(s)
        t = getfirst(x->occursin("ssh -f $(ssh_string(s))", x), split(read(pipeline(`ps aux` , stdout = `grep Resource`), String), "\n"))
        if t !== nothing
            run(`kill $(split(t)[2])`)
        end
    end
        
    if s.local_port != 0
        destroy_tunnel(s)
    end

    # Here we check what the modify time of the server-side localhost file is.
    # The server will rewrite the file with the correct port, which we use to see
    # whether the server started succesfully.
    function checktime()
        curtime = 0
        try
            cmd = "stat -c %Z  ~/.julia/config/DFControl/servers/localhost.json"
            curtime = islocal(s) ? mtime(config_path("servers", "localhost.json")) : parse(Int, server_command(s.username, s.domain, cmd)[1])
        catch
            nothing
        end
        return curtime
    end
    firstime = checktime()

    # Now we upload the user's uuid for safety
    if !islocal(s)
        push(config_path("user_uuid"), s, "~/.julia/config/DFControl/user_uuid")
    end

    julia_cmd = """$(s.julia_exec) --startup-file=no -t auto -e "using DFControl; DFControl.Resource.run()" &> ~/.julia/config/DFControl/logs/daemon.log"""
    if s.domain != "localhost"
        run(Cmd(`ssh -f $(ssh_string(s)) $julia_cmd`, detach=true))
    else
        scrpt = "using DFControl; DFControl.Resource.run()"
        e = s.julia_exec
        julia_cmd = `$(e) --startup-file=no -t auto -e $(scrpt) '&''>' '~'/.julia/config/DFControl/errors.log '&'`
        run(Cmd(julia_cmd, detach=true), wait=false)
    end
        
    #TODO: little hack here
    retries = 0
    while checktime() <= firstime && retries < 60
        retries += 1
        sleep(1)
    end

    if retries == 60
        error("Something went wrong starting the server.")
    else
        if s.local_port == 0
            @info "Daemon on Server $(s.name) started, listening on port $(s.port)."
        else
            tserver = islocal(s) ? read_config(config_path("servers/localhost.json")) : load_remote_config(s.username, s.domain)
            s.port = tserver.port
            construct_tunnel(s)
            @info "Daemon on Server $(s.name) started, listening on local port $(s.local_port)."
            @info "Saving updated server info..."
            save(s)
        end
    end
     
    return s
end

function maybe_start(s::Server)
    if !isalive(s)
        error("Server $(s.name) is not running.\nTry to start it with Servers.start(Server($(s.name)))")
    end
    return s
end
maybe_start(j) = (s = Server(j); return maybe_start(s))

"""
    kill(s::Server)

Kills the daemon process on [`Server`](@ref) `s`.
"""
kill(s::Server) = HTTP.put(s, "/kill_server")

function restart(s::Server)
    kill(s)
    return start(s)
end

function ask_input(::Type{T}, message, default=nothing) where {T}
    if default === nothing
        t = ""
        print(message * ": ")
        while isempty(t)
            t = readline()
        end
    else
        print(message * " (default: $default): ")
        t = readline()
        if isempty(t)
            return default
        end
    end
    if T != String
        return parse(T, t)
    else
        return t
    end
end
   
function maybe_create_localhost()
    t = load("localhost")
    if t === nothing
        @info "Initializing localhost server configuration."
        scheduler = Sys.which("sbatch") === nothing ? Bash : Slurm
        if !haskey(ENV, "DFCONTROL_PORT")
            port = ask_input(Int, "Port", 8080)
        else
            port = parse(Int, ENV["DFCONTROL_PORT"])
        end
        if !haskey(ENV, "DFCONTROL_JOBDIR")
            dir = ask_input(String, "Default Jobs directory", homedir())
            if dir != homedir()
                while !ispath(dir)
                    @warn "$dir, no such directory"
                    dir = ask_input(String, "Default Jobs directory", homedir())
                end
            end
        else
            dir = ENV["DFCONTROL_JOBDIR"]
        end
        julia_exec = joinpath(Sys.BINDIR, "julia")
        max_concurrent_jobs = ask_input(Int, "Max Concurrent Jobs", 100)
        out = Server("localhost", ENV["USER"], "localhost", port, scheduler, "", julia_exec,
                     dir, 0, max_concurrent_jobs)
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
    if islocal(server)
        cp(server_file, filename; force = true)
    else
        out = Pipe()
        err = Pipe()
        run(pipeline(`scp $(ssh_string(server) * ":" * server_file) $filename`, stdout=out, stderr=err))
        close(out.in)
        close(err.in)
    end
end

function pull(j::Job, f, t)
    @assert ispath(j, f) "File $f not found in jobdir."
    pull(Server(j.server), joinpath(j, f), t)
end

"""
    push(local_file::String, server::Server, server_file::String)

Pushes the `local_file` to the `server_file` on the server.
"""
function push(filename::String, server::Server, server_file::String)
    if islocal(server)
        cp(filename, server_file; force = true)
    else
        run(`scp $filename $(ssh_string(server) * ":" * server_file)`)
    end
end

"Executes a command through `ssh`."
function server_command(username, domain, cmd)
    out = Pipe()
    err = Pipe()
    if domain == "localhost"
        process = run(pipeline(ignorestatus(cmd), stdout=out, stderr=err))
    else
        process = run(pipeline(ignorestatus(`ssh $(username * "@" * domain) $cmd`), stdout=out, stderr=err))
    end
    close(out.in)
    close(err.in)

    stdout = read(out, String)
    stderr = read(err, String)
    return (
      stdout = stdout,
      stderr = stderr,
      exitcode = process.exitcode
    )
end
    
server_command(s::Server, cmd) = server_command(s.username, s.domain, cmd)

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

function Base.readdir(s::Server, dir::String)
    maybe_start(s)
    resp = HTTP.get(s, "/readdir/" * abspath(s, dir))
    return JSON3.read(resp.body, Vector{String})
end

Base.abspath(s::Server, p) =
    isabspath(p) ? p : joinpath(s, p)

function Base.mtime(s::Server, p)
    resp = HTTP.get(s, "/mtime/" * p)
    JSON3.read(resp.body, Float64)
end
end
