using DFControl, Test
DFControl.Service.configure_pseudoset("test", joinpath(testdir, "testassets", "pseudos"))
@test DFControl.Service.pseudos("test", "UPF")[:Ni] ==
      Pseudo("Ni.pbesol-n-kjpaw_psl.0.1.UPF", joinpath(testdir, "testassets", "pseudos"), 41., 236.)
