using DFControl, Test

Client.rm_pseudoset!(test_server, "test")

@test_throws DFC.Client.HTTP.ExceptionRequest.StatusError Client.pseudos(test_server, "test", [:Ni])

rm(test_server, Environment("test_default"))
