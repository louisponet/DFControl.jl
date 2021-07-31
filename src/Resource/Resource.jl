module Resource
    using HTTP, JSON3
    using ..DFControl, ..Service
    
    const ROUTER = HTTP.Router()

    save_job(req) = Service.save_job(JSON3.read(req.body, DFJob))
    HTTP.@register(ROUTER, "POST", "/jobs", save_job)

    get_job(req) = Service.get_job(joinpath(HTTP.URIs.splitpath(req.target)[2:end]...))
    HTTP.@register(ROUTER, "GET", "/jobs/*", get_job)

    kill_server(req) = exit()
    HTTP.@register(ROUTER, "PUT", "/kill_server", kill_server)

    function requestHandler(req)
        obj = HTTP.handle(ROUTER, req)
        return HTTP.Response(200, JSON3.write(obj))
    end

    function run(port)
        @async HTTP.serve(requestHandler, "0.0.0.0", port)
        # Service.start()
        Service.main_loop()
    end    
end
