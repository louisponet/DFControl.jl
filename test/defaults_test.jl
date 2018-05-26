
using DFControl

prevdefault = getdefault_server()
setdefault_server("localhost")
setdefault_pseudodir(:test, joinpath(@__DIR__, "testassets/pseudos"))
# setdefault_pseudodir(:test, joinpath(Pkg.dir("DFControl"), "test/testassets/pseudos"))
configuredefault_pseudos(pseudo_dirs=Dict(:test => getdefault_pseudodirs()[:test]))

@add_default testdefaultstring = "test"

@test testdefaultstring == "test"
@test DFControl.getdefault_pseudo(:Pt, :test) == "Pt.UPF"
@test DFControl.getdefault_server() == "localhost"
