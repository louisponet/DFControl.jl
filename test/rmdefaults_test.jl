using DFControl, Test

DFControl.removedefault_pseudos(:test)

@test DFControl.getdefault_pseudo(:Si, :test) == nothing

removedefault(:testdefaultstring)
@test DFControl.testdefaultstring == nothing
setdefault_server(prevdefault)
