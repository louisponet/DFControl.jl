sshcmd(server, cmd) = run(`ssh -t $server $cmd`)
sshreadstring(server, cmd) = read(`ssh -t $server $cmd`, String)
"Tests whether a directory exists on a server and if not, creates it."
function mkserverdir(server, dir)
    testfile = joinpath(dir, "tmp.in")
    try
        run(`ssh -t $server touch $testfile`)
        run(`ssh -t $server rm $testfile`)
    catch
        run(`ssh -t $server mkdir -p $dir`)
        @info "$dir did not exist on $server, it was created."
    end
end

Base.joinpath(s::Server, p...) = joinpath(s.default_jobdir, p...)
#Gives the reverse (last job is listed first) of the output, omitting the header lines

