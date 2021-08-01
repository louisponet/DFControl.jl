"""
    qstat(server)

If sbatch is running on server it will return the output of the `qstat` command.
"""
qstat(server) = server == "localhost" ? run(`qstat`) : sshcmd(server, "qstat")
qstat() = qstat(getdefault_server())

"""
    watch_qstat(server)
If sbatch is running on server it will return the output of the `watch qstat` command.
"""
function watch_qstat(server)
    return server == "localhost" ? run(`watch qstat`) : sshcmd(server, "watch qstat")
end
watch_qstat() = watch_qstat(getdefault_server())

#--------------- Slurm Interactions -------------------#
