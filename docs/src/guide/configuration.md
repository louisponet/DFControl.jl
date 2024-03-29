# Configuration

## [Servers](@ref servers_header)
Since `DFControl` utilizes a client-server rest-api model, each server will have its own local deamon running, which stores
certain server-side items such as pseudopotentials, [`Environments`](@ref Environment) and [`Execs`](@ref Exec).

Make sure that both `julia` and `DFControl` are installed on the target remote server before trying to create a connection to it.
Also, it is necessary that your ssh keys are registered on the remote so that no passwords are required for `ssh` connections.

To up a new [`Server`](@ref), in this case with `name = "daint"` an interactive menu can be called like:
```julia
Server("daint")
```

This will set up all the required information and store it for later use. To start a [`Server`](@ref) simply call `start(server)`, which
will launch the daemon externally.
Make sure that the `host` address points to a persistent server, i.e. if a cluster has multiple frontend nodes, find the ip address of one particular one,
and set that one as the `host` address.
During the setup of an external server, it may ask whether a `local_tunnel` should be created. If enabled, a persistent ssh tunnel will be
created from the local machine to the target host through which the http requests will be sent. This is useful when the remote is behind an
authentication firewall.

To reconfigure a [`Server`](@ref) you can always do:
```julia
s = Server("daint")
s.root_jobdir = "<where you want jobs to be saved>"
save(s)
```
If a [`Server`](@ref) was running first stop it using `kill(server)`.

## Pseudopotentials
Pseudopotentials are grouped in sets, which are stored for later ease of use.
They can be set up using the [`configure_pseudoset`](@ref) function.
```julia
configure_pseudoset(local_server(), "set_name", "/dir/to/pseudos")
```

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
Here we will set up an environment and save it on the local [`Server`](@ref). Change the information according to your own setup.
```julia
e = Environment(name="default", parallel_exec=Exec(exec="mpirun", flags=Dict("-np"=> 4), directives=Dict("N" => 1, "partition" => "parallel"), exports=Dict("OMP_NUM_THREADS"=>1)))
save(local_server(), e)
```

## [Execs](@ref execs_header)
An [`Exec`](@ref) embodies an executable as the name suggests.
They hold both the executable, directory and modules required to run the executable.
Similar to [`Environments`](@ref Environment), they are stored on the server for later use. When storing it is verified that they are able to run.

```julia
e = Exec(name="pw", exec="pw.x", dir="/home/user/Softare/qe/bin", modules = ["intel", "intel-mpi", "intel-mkl"])
save(local_server(), e)
```


## Loading/Saving
As shown above there are three main methods related to storing and retrieving data i.e.:
```@docs
DFControl.Client.RemoteHPC.load
DFControl.Client.RemoteHPC.save
```

After completing this configuration, it is suggested to look at the (Basic Usage)[@ref] for an introduction to the basic usage of `DFControl`.
