using DFControl, Test

configure_pseudoset(test_server,"test", joinpath(testdir, "testassets", "pseudos"))

save(test_server, DFControl.Client.environment_from_jobscript(test_server, "test_default", joinpath(testdir, "testassets", "reference_job", "job.tt")))
