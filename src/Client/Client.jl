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
using ..DFControl: config_path

using ..Calculations: set_name!, set_kpoints!, data; export set_name!, set_kpoints!, data
using ..Structures: set_pseudos!, element; export set_pseudos!, element
using ..Jobs: set_flow!; export set_flow!
import ..DFControl: bandgap

include("job.jl")
export submit, save, isrunning, state, versions, last_version, switch_version!, rm_version!, abort,
       outputdata, readfermi, readbands, bandgap, archive, cleanup
include("pseudos.jl")
export configure_pseudoset, rm_pseudoset!, list_pseudosets
# include("firecrest.jl")

end
