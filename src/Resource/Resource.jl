module Resource
    using HTTP, JSON3
    using ..DFControl, ..Service
    using ..Utils
    
    const ROUTER = HTTP.Router()

    include("job.jl")
    
    kill_server(req) = exit()
    HTTP.@register(ROUTER, "PUT", "/kill_server", kill_server)

    get_server_config(req) = Service.server_config()
    HTTP.@register(ROUTER, "GET", "/server_config", get_server_config)

    function requestHandler(req)
        obj = HTTP.handle(ROUTER, req)
        if obj === nothing
            return HTTP.Response(204)
        else
            return HTTP.Response(200, JSON3.write(obj))
        end
    end

    function run(port)
        @async HTTP.serve(requestHandler, "0.0.0.0", port)
        # Service.start()
        Service.main_loop()
    end    
end
