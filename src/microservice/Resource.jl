module Resource
    using HTTP, JSON3, StructTypes
    using ..DFControl, ..Service
    @inline function JSON3.read(::StructTypes.Mutable, buf, pos, len, b, ::Type{DFCalculation}; kw...)
        x = DFCalculation{DFControl.NoPackage}("", execs=Exec[])
        pos, x = JSON3.read!(StructTypes.Mutable(), buf, pos, len, b, DFCalculation, x; kw...)
        return pos, x
    end

    const ROUTER = HTTP.Router()

    save_job(req) = Service.save_job(JSON3.read(req.body, DFJob))
    HTTP.@register(ROUTER, "POST", "/jobs", save_job)

    get_job(req) = Service.get_job(joinpath(HTTP.URIs.splitpath(req.target)[2:end]...))
    HTTP.@register(ROUTER, "GET", "/jobs/*", get_job)

    function requestHandler(req)
        obj = HTTP.handle(ROUTER, req)
        return HTTP.Response(200, JSON3.write(obj))
    end

    function run()
        HTTP.serve(requestHandler, "0.0.0.0", 8080)
    end    
end
