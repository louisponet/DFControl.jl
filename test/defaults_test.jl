
using DFControl

setdefault_server("localhost")
setdefault_pseudodir(:test, joinpath(Pkg.dir("DFControl"), "test/testassets/pseudos"))
configure_defaultpseudos()


@add_default testdefaultstring = "test"

@test testdefaultstring == "test"
@test DFControl.getdefault_pseudo(:Si, :test) == "Si.UPF"
@test DFControl.getdefault_server() == "localhost"
