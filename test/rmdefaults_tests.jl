using DFControl, Test

Client.rm_pseudoset!("test", "localhost_test")

@test_throws ErrorException Client.pseudos("localhost_test", "test", [:Ni])
