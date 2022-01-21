config_path(args...) = joinpath(abspath(first(DEPOT_PATH), "config", "DFControl"), args...)

# All the folder
mkpath(config_path())
mkpath(config_path("servers"))
mkpath(config_path("execs"))
mkpath(config_path("jobs"))
mkpath(config_path("logs"))
mkpath(config_path("logs/jobs"))

# All the job info
touch(config_path("jobs", "pending.txt"))
touch(config_path("jobs", "archived.txt"))
touch(config_path("jobs", "active.txt"))
touch(config_path("jobs", "running.txt"))
touch(config_path("workflows", "pending.txt"))
touch(config_path("workflows", "archived.txt"))
touch(config_path("workflows", "active.txt"))
touch(config_path("workflows", "running.txt"))
