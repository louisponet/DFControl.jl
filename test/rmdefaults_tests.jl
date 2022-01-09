using DFControl, Test

Client.rm_pseudoset!("test", server= "localhost")

@test_throws ErrorException Client.pseudos("localhost", "test", [:Ni])
