# Configuration

## [Servers](@ref servers_header)
Since `DFControl` utilizes a client-server rest-api model, each server will have its own local deamon running, which stores
certain server-side items such as pseudopotentials, [`Environments`](@ref Environment) and [`Execs`](@ref Exec).
To set up a [`Server`](@ref), an interactive menu can be called like:
```julia
Server("daint")
```

This will set up all the required information and store it for later use. To start a [`Server`](@ref) simply call `start(server)`, which
will launch the daemon externally.


## Pseudopotentials
Pseudopotentials are grouped in sets, which are stored for later ease of use.
They can be set up using the [`configure_pseudoset`](@ref) function.
This will go through the specified directories and find files that follow the naming convention
`element.x_y_z.xyz` e.g. `Si.pbesol-n-kjpaw_psl.0.1.UPF` or `si.pbesol-n-kjpaw_psl.0.1.UPF`. 
If multiple are found, all will be stored and the required one can be later specified.

The pseudos will remain on the server where they are stored, and can be listed using [`list_pseudosets`](@ref).
See [Pseudo Potentials](@ref pseudo_header) for further usage details.

```@docs
configure_pseudoset
list_pseudosets
rm_pseudoset!
```

## [Environments](@ref environments_header)
Environments specify the skeleton of the job script, i.e. which environment variables need to be set, which scheduler flags, etc.
Here we will set up an environment and save it on server `daint`.
```julia
e = Environment(name="normal_1node", MPI_command="srun", scheduler_flags=["#SBATCH -N 1", "#SBATCH --partition=parallel"], exports=["OMP_NUM_THREADS=1"])
save(Server("daint"), e)
```
Alternatively, one can extract an environment from a preexisting jobscript using
```@docs
DFC.Client.environment_from_jobscript
```

## [Execs](@ref execs_header)
An [`Exec`](@ref) embodies an executable as the name suggests.
They hold both the executable, directory and modules required to run the executable.
Similar to [`Environments`](@ref Environment), they are stored on the server for later use. When storing it is verified that they are able to run.

```julia
e = Exec(name="pw", exec="pw.x", dir="/home/user/Softare/qe/bin", modules = ["intel", "intel-mpi", "intel-mkl"])
save(Server("daint"), e)
```


## Loading/Saving
As shown above there are three main methods related to storing and retrieving data i.e.:
```@docs
Database.load
Database.save
rm(::Database.Storable)
```

After completing this configuration, it is suggested to look at the (Basic Usage)[@ref] for an introduction to the basic usage of `DFControl`.
