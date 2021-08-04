module Client
    using HTTP, JSON3, StructTypes, Dates, JLD2, Distributed, REPL.TerminalMenus
    using ..DFControl
    using ..Utils
    using ..Servers
    using ..Calculations
    using ..Jobs
    
    @inline function JSON3.read(::StructTypes.Mutable, buf, pos, len, b, ::Type{Calculation}; kw...)
        x = Calculation{DFControl.NoPackage}("", execs=Exec[])
        pos, x = JSON3.read!(StructTypes.Mutable(), buf, pos, len, b, Calculation, x; kw...)
        p = eltype(x.execs)
        t = Calculation{p}(x.name, x.dir, x.flags, x.data, x.execs, x.run, x.outdata, x.infile, x.outfile)
        return pos, t
    end
   
   

    # include("job.jl")
    # export save, isrunning, versions, last_version
    
    # include("pseudos.jl")
    
        
end
