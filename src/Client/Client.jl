module Client
    using HTTP, JSON3, StructTypes, Dates, JLD2, Distributed, REPL.TerminalMenus
    using ..DFControl
    using ..DFControl: Server
    using ..Utils
    
    function package(execs::Vector{Exec})
        if any(DFC.is_qe_exec, execs)
            return QE
        elseif any(DFC.is_wannier_exec, execs)
            return Wannier90
        elseif any(DFC.is_elk_exec, execs)
            return Elk
        elseif any(DFC.is_abinit_exec, execs)
            return Abinit
        end
    end

    @inline function JSON3.read(::StructTypes.Mutable, buf, pos, len, b, ::Type{Calculation}; kw...)
        x = Calculation{DFControl.NoPackage}("", execs=Exec[])
        pos, x = JSON3.read!(StructTypes.Mutable(), buf, pos, len, b, Calculation, x; kw...)
        p = eltype(x.execs)
        t = Calculation{p}(x.name, x.dir, x.flags, x.data, x.execs, x.run, x.outdata, x.infile, x.outfile)
        return pos, t
    end
   
    HTTP.request(method::String, s::Server, url, args...; kwargs...) =
        HTTP.request(method, string(http_string(s), url), args...; kwargs...)
        
    for f in (:get, :put, :post, :head)
        str = uppercase(string(f))
        @eval HTTP.$(f)(s::Server, url, args...; kwargs...) =
            HTTP.request("$($str)", s, url, args...; kwargs...)
    end


    include("job.jl")
    include("calculation.jl")
    export save, isrunning, versions, last_version
    
    include("server.jl")
    include("pseudos.jl")
    
        
end
