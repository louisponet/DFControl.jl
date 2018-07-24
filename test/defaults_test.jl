
using DFControl
import DFControl: Exec, data, add!, execs

prevdefault = getdefault_server()
setdefault_server("localhost")
testdir = joinpath(Pkg.dir("DFControl"), "test/")
setdefault_pseudodir(:test, joinpath(testdir, "testassets/pseudos/"))
# setdefault_pseudodir(:test, joinpath(Pkg.dir("DFControl"), "test/testassets/pseudos"))
configuredefault_pseudos(pseudo_dirs=Dict(:test => getdefault_pseudodirs()[:test]))
default_pseudos[:Pt]
@add_default testdefaultstring = "test"

@test testdefaultstring == "test"
@test DFControl.getdefault_pseudo(:Pt, :test) == "Pt.UPF"
@test DFControl.getdefault_server() == "localhost"
