using DFControl, Test
configure_pseudoset("test", joinpath(testdir, "testassets", "pseudos"), server="localhost")
add_environment(DFControl.Client.environment_from_jobscript(joinpath(testdir, "testassets", "reference_job", "job.tt"), server="localhost"), "test_default", server="localhost")
# @test DFControl.Service.pseudos("test", "UPF")[:Ni] ==
#       Pseudo("Ni.pbesol-n-kjpaw_psl.0.1.UPF", joinpath(testdir, "testassets", "pseudos"), 41., 236.)
