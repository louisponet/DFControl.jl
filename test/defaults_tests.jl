using DFControl, Test
configure_pseudoset("test", joinpath(testdir, "testassets", "pseudos"), "localhost_test")
add_environment(DFControl.Client.environment_from_jobscript(joinpath(testdir, "testassets", "reference_job", "job.tt"), "localhost_test"), "test_default", "localhost_test")
# @test DFControl.Service.pseudos("test", "UPF")[:Ni] ==
#       Pseudo("Ni.pbesol-n-kjpaw_psl.0.1.UPF", joinpath(testdir, "testassets", "pseudos"), 41., 236.)
