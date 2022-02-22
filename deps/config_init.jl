using UUIDs
config_path(args...) = joinpath(abspath(first(DEPOT_PATH), "config", "DFControl"), args...)

# All the folder
mkpath(config_path())
mkpath(config_path("jobs"))
mkpath(config_path("workflows"))
mkpath(config_path("logs"))
mkpath(config_path("logs/daemon/"))
mkpath(config_path("logs/jobs"))
mkpath(config_path("logs/servers"))
mkpath(config_path("logs/jobs/$(gethostname())"))
mkpath(config_path("logs/servers/$(gethostname())"))
mkpath(config_path("storage"))
mkpath(config_path("storage/servers"))
mkpath(config_path("storage/execs"))
mkpath(config_path("storage/pseudos"))
mkpath(config_path("storage/environments"))

# All the job info
touch(config_path("jobs", "pending.txt"))
touch(config_path("jobs", "archived.txt"))
touch(config_path("jobs", "active.txt"))
touch(config_path("jobs", "running.txt"))
touch(config_path("workflows", "pending.txt"))
touch(config_path("workflows", "archived.txt"))
touch(config_path("workflows", "active.txt"))
touch(config_path("workflows", "running.txt"))
