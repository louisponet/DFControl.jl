module Client
using ..Database
import ..Database: load, save
using HTTP, JSON3, StructTypes, Dates, JLD2, Distributed, REPL.TerminalMenus, Reexport
using ..Utils
using ..FileIO

@reexport using ..Servers
@reexport using ..Structures
@reexport using ..Calculations
@reexport using ..Jobs

# @inline function JSON3.read(::StructTypes.Mutable, buf, pos, len, b, ::Type{Calculation};
#                             kw...)
#     x = Calculation{Calculations.NoPackage}(""; exec = Exec())
#     pos, x = JSON3.read!(StructTypes.Mutable(), buf, pos, len, b, Calculation, x; kw...)
#     if x.exec.exec ∈ Calculations.WAN_EXECS
#         p = Wannier90
#     elseif x.exec.exec ∈ Calculations.QE_EXECS
#         p = QE
#     elseif x.exec.exec ∈ Calculations.ELK_EXECS
#         p = ELK
#     else
#         @warn "Package not identified from execs $(x.exec)."
#     end

#     for d in x.data
#         if d.data[1] isa Vector
#             dat = NTuple{length(d.data[1]), Float64}[]
#             for d_ in d.data
#                 push!(dat, (d_...,))
#             end
#             d.data = dat
#         end
#     end

#     t = Calculation{p}(x.name, x.flags, x.data, x.exec, x.run,  x.infile,
#                        x.outfile)
#     return pos, t
# end

using ..Calculations: set_name!, set_kpoints!, data; export set_name!, set_kpoints!, data
using ..Structures: set_pseudos!, element; export set_pseudos!, element
using ..Jobs: set_flow!; export set_flow!
import ..DFControl: bandgap

include("job.jl")
export submit, save, isrunning, state, versions, last_version, switch_version!, rm_version!, abort,
       outputdata, readfermi, readbands, bandgap, archive, cleanup
include("pseudos.jl")
export configure_pseudoset, rm_pseudoset!, list_pseudosets

end
