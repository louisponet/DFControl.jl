
using DFControl, Test
import DFControl: Exec, data, add!, execs

prevdefault = getdefault_server()
setdefault_server("localhost")
setdefault_pseudodir(:test, joinpath(testdir, "testassets", "pseudos"))
configuredefault_pseudos(pseudo_dirs=Dict(:test => getdefault_pseudodirs()[:test]))

@test DFControl.getdefault_pseudo(:Pt, :test) == "Pt.UPF"
@test DFControl.getdefault_server() == "localhost"
