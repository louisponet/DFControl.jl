const SERVER_DIR = DFC.config_path("servers")
const LOCALHOST_UP = Ref(false)

function known_servers(fuzzy="")
    if ispath(SERVER_DIR)
        servers = [JLD2.load(joinpath(SERVER_DIR, s))["server"] for s in filter(x -> occursin(fuzzy, x), readdir(SERVER_DIR))]
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

function establish_connection(server::Server)
    if server.name == "localhost"
        proc = Distributed.addprocs(1)[1]
    else
        proc = Distributed.addprocs([(ssh_string(server), 1)])[1]
    end
    return proc
end

remove_server!(name::String) = ispath(joinpath(SERVER_DIR, name * ".jld2")) && rm(joinpath(SERVER_DIR, name * ".jld2"))

function start(s::Server)
    @info "Starting:\n$s"
    cmd = Cmd(`$(s.julia_exec) --startup-file=no -t auto -e "using DFControl; DFControl.Resource.run($(s.port))"`; detach = true)
    proc = establish_connection(s)
    
    p   = remotecall(run, proc, cmd; wait = false)
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
    rmprocs(p)
end

function maybe_start_server(s::Server)
    try
        conf = JSON3.read(HTTP.get(s, "/server_config").body, Server)
        # @assert s.scheduler == conf.scheduler "Scheduler mismatch."
    catch
        start(s)
    end
    return s
end
maybe_start_server(j::Union{Job, String}) = (s = Server(j); return maybe_start_server(s))

kill_server(s::Server) = HTTP.put(s, "/kill_server")

function restart_server(s::Server)
    kill_server(s)
    start(s)
end

function maybe_create_localhost()
    t = load_server("localhost")
    if t === nothing
        scheduler = Sys.which("sbatch") === nothing ? DFC.Bash : DFC.Slurm
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
        out = Server("localhost", ENV["USER"], "localhost", port, scheduler, "", julia_exec, dir)
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
        cp(server_file, filename, force=true)
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
        cp(filename, server_file, force=true)
    else
        run(`scp $filename $(ssh_string(server) * ":" * server_file)`)
    end
end


