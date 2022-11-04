using DFControl, Test

configure_pseudoset(test_server,"test", joinpath(testdir, "testassets", "pseudos"))


e = Environment("test", Dict("N" => 1, "time" => "00:01:01"),
                            Dict("OMP_NUM_THREADS" => 1), "", "",
                            Exec(; name = "srun", exec = "srun"))
save(test_server, e)
save(test_server, Exec(; name = "pw", exec = "pw.x"))
