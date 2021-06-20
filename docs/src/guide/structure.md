# Structure
## Contents
```@contents
Pages=["structure.md"]
```
## Index
```@index
Pages=["structure.md"]
```


```@docs
Structure
Base.getindex(::Structure, ::Int)
update_geometry!
polyhedron
```

## Symmetries
This functionality relies on `spglib` to find the symmetries of the `Structure` and 
supply various related quantities.
```@docs
high_symmetry_kpath
high_symmetry_kpoints
international
niggli_reduce
symmetry_operators
```

## Cell
!!! note
    The lattice vectors are stored as the *columns* of the `cell`.
    
```@docs
a
b
c
cell_parameters
volume
create_supercell
scale_cell!
```

## Atom
```@docs
Atom
atoms(::Structure)
set_position!
```
## Element
```@docs
DFControl.Element
element
```

## [Pseudo Potentials](@id pseudo_header)
```@docs
Pseudo
set_pseudo!
set_pseudos!
```
!!! note
If pseudopotentials from different sets/directories are specified for the atoms, they will be 
copied to the `DFJob` `local_dir`.

### Pseudo sets
See the section on [Configuration](@ref) for a demonstration on how to set up pseudopotential sets.
```@docs
getdefault_pseudo
setdefault_pseudodir
configuredefault_pseudos
```

## Magnetization
Magnetization can be set on a per [`Atom`](@ref Atom) basis. It will partially determine the unique
[`Atom`](@ref Atom) types and also what calculation flags should be set in order to allow for the 
magnetic calculation (either colinear or non-colinear).
```@docs
set_magnetization!
```

## Projections
These projections will mainly be used for generating Wannier90 inputs, and to distinguish which indices
in the Wannier90 output matrices correspond to the various atoms/orbitals.
```@docs
Projection
set_projections!
```

### Orbitals
```@docs
Orbital
orbital
```

## DFT + U 
```@docs
DFTU
set_Hubbard_U!
set_Hubbard_J0!
set_Hubbard_α!
set_Hubbard_β!
set_Hubbard_J!
```
