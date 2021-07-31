module Client
    using HTTP, JSON3, StructTypes
    using ..DFControl
    using ..DFControl: Server

    @inline function JSON3.read(::StructTypes.Mutable, buf, pos, len, b, ::Type{DFCalculation}; kw...)
        x = DFCalculation{DFControl.NoPackage}("", execs=Exec[])
        pos, x = JSON3.read!(StructTypes.Mutable(), buf, pos, len, b, DFCalculation, x; kw...)
        return pos, x
    end
    
    http_string(s::Server) = "http://$(s.domain):$(s.port)"

    HTTP.request(method::String, s::Server, url, args...) =
        HTTP.request(method, string(http_string(s), url), args...)
        
    for f in (:get, :put, :post, :head)
        @eval HTTP.$(f)(s::Server, url, args...) =
            HTTP.request(uppercase(string($f)), s, url, args...)
    end

    function DFControl.DFJob(s::Server, dir::String)
        job = JSON3.read(HTTP.get(s, joinpath("/jobs", dir)).body, DFJob)
        job.server = s.name
        job.server_dir = job.local_dir
        return job
    end

    kill_server(s::Server) = HTTP.put(s, "/kill_server")

        
        
end
