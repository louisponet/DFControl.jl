module Client
using HTTP, JSON3, StructTypes, Dates, JLD2, Distributed, REPL.TerminalMenus, Reexport
using ..DFControl
using ..Utils
using ..FileIO
@reexport using ..Servers
@reexport using ..Structures
@reexport using ..Calculations
@reexport using ..Jobs

@inline function JSON3.read(::StructTypes.Mutable, buf, pos, len, b, ::Type{Calculation};
                            kw...)
    x = Calculation{Calculations.NoPackage}(""; execs = Exec[])
    pos, x = JSON3.read!(StructTypes.Mutable(), buf, pos, len, b, Calculation, x; kw...)
    if any(y -> y.exec ∈ Calculations.WAN_EXECS, x.execs)
        p = Wannier90
    elseif any(y -> y.exec ∈ Calculations.QE_EXECS, x.execs)
        p = QE
    elseif any(y -> y.exec ∈ Calculations.ELK_EXECS, x.execs)
        p = ELK
    else
        @warn "Package not identified from execs $(x.execs)."
    end

    t = Calculation{p}(x.name, x.dir, x.flags, x.data, x.execs, x.run, x.outdata, x.infile,
                       x.outfile)
    return pos, t
end

using ..DFControl: set_dir!; export set_dir!
using ..Calculations: set_name!, set_kpoints!, data; export set_name!, set_kpoints!, data
using ..Structures: set_pseudos!, element; export set_pseudos!, element


include("job.jl")
export submit, save, isrunning, versions, last_version, outputdata

include("pseudos.jl")
export configure_pseudoset, rm_pseudoset!

end
