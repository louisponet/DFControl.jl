# Configuration

In order to not always have to respecify where to look for PseudoPotential sets, 
they can be set up once and be remembered by DFControl. 
This can be done through:

```julia
using DFControl

setdefault_pseudodir(:pseudo_set_name1, "/absolute/path/to/pseudopotential/set/1")
setdefault_pseudodir(:pseudo_set_name2, "/absolute/path/to/pseudopotential/set/2")

configuredefault_pseudos()
```
This will go through the specified directories and find files that follow the naming convention
`element.x_y_z.xyz` e.g. `Si.pbesol-n-kjpaw_psl.0.1.UPF` or `si.pbesol-n-kjpaw_psl.0.1.UPF`. 
If multiple are found, all will be stored and the required one can be later specified.
See [Pseudo Potentials](@ref pseudo_header) for further usage details.

In order to change the pseudos associated with a set name simply:
```julia
setdefault_pseudodir(:pseudo_set_name1, "/new/absolute/path/to/pseudopotential/set/1")
setdefault_pseudodir(:pseudo_set_name2, "/new/absolute/path/to/pseudopotential/set/2")

configuredefault_pseudos()
```

In order to remove a default pseudo set one can run:
```julia
removedefault_pseudos(:pseudo_set_name)
```

After completing this configuration, it is suggested to look at the (Basic Usage)[@ref] for an introduction to the basic usage of DFControl.
