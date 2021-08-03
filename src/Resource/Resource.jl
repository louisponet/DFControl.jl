module Resource
    # This module handles all the handling, parsing and transforming HTTP commands,
    # and brokers between HTTP requests and Service that fullfils them. 
    using HTTP, JSON3, Dates
    using ..DFControl, ..Service
    using ..Utils
    
    const ROUTER = HTTP.Router()

    include("job.jl")

    # GENERAL
    kill_server(req) = exit() 
    HTTP.@register(ROUTER, "PUT", "/kill_server", kill_server)

    get_server_config(req) = Service.server_config()
    HTTP.@register(ROUTER, "GET", "/server_config", get_server_config)

    get_ispath(req) = ispath(job_path(req))
    HTTP.@register(ROUTER, "GET", "/get_ispath/*", get_ispath)
    
    # PSEUDOS
    
    function pseudos(req)
        fuzzy = String(req.body)
        fuzzy = fuzzy == "\"\"" ? "" : fuzzy
        return Service.pseudos( splitpath(req.target)[end], fuzzy)
    end
    HTTP.@register(ROUTER, "GET", "/pseudos/*", pseudos)
    
    pseudo_sets(req) = Service.pseudo_sets()
    HTTP.@register(ROUTER, "GET", "/pseudo_sets/", pseudo_sets)
    
    configure_pseudos(req) = Service.configure_pseudos(req.body, job_path(req))
    HTTP.@register(ROUTER, "POST", "/configure_pseudos/*", configure_pseudos)
    
    rm_pseudos!(req) = Service.rm_pseudos!(req.body)
    HTTP.@register(ROUTER, "PUT", "/rm_pseudos", rm_pseudos!)
    # EXECS
    verify_exec(req) = Service.verify_exec(JSON3.read(req.body, Exec))
    HTTP.@register(ROUTER, "GET", "/verify_exec", verify_exec)
    
    # RUNNING
    
    function requestHandler(req)
        start = Dates.now()
        @info (timestamp=start, event="ServiceRequestBegin", tid=Threads.threadid(), method=req.method, target=req.target, body=req.body)
        local resp
        try
            # DFC.Revise.revise()
            # obj = Base.invokelatest(HTTP.handle,ROUTER, req)
            obj = HTTP.handle(ROUTER, req)
            if obj === nothing
                resp = HTTP.Response(204)
            else
                resp = HTTP.Response(200, JSON3.write(obj))
            end
        catch e
            s = IOBuffer()
            showerror(s, e, catch_backtrace(); backtrace=true)
            errormsg = String(resize!(s.data, s.size))
            @error errormsg
            resp = HTTP.Response(500, errormsg)
        end
        stop = Dates.now()
        @info (timestamp=stop, event="ServiceRequestEnd", tid=Threads.threadid(), method=req.method, target=req.target, duration=Dates.value(stop - start), status=resp.status, bodysize=length(resp.body))
        return resp
    end

    function run(port)
        cd(DFC.Server("localhost").default_jobdir)
        Service.global_logger(Service.daemon_logger()) 
        # Service.start()
        # server = HTTP.Sockets.listen(HTTP.Sockets.InetAddr(parse(IPAddr, "0.0.0.0"), port))
        @async HTTP.serve(requestHandler, "0.0.0.0", port)
        Service.main_loop()
    end    
end