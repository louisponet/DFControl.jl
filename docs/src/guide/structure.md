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
Structures.update_geometry!
Structures.polyhedron
```

## Symmetries
This functionality relies on `spglib` to find the symmetries of the `Structure` and 
supply various related quantities.
```@docs
Structures.high_symmetry_kpath
Structures.high_symmetry_kpoints
Structures.international
Structures.niggli_reduce
Structures.symmetry_operators
```

## Cell
!!! note
    The lattice vectors are stored as the *columns* of the `cell`.
    
```@docs
Structures.a
Structures.b
Structures.c
Structures.cell_parameters
Structures.volume
Structures.create_supercell
Structures.scale_cell!
```

## Atom
```@docs
Atom
Structures.set_position!
```
## Element
```@docs
Structures.Element
element
```

## [Pseudo Potentials](@id pseudo_header)
```@docs
set_pseudos!
configure_pseudoset
list_pseudosets
```
See the section on [Configuration](@ref) for a demonstration on how to set up pseudopotential sets.

## Magnetization
Magnetization can be set on a per [`Atom`](@ref Atom) basis. It will partially determine the unique
[`Atom`](@ref Atom) types and also what calculation flags should be set in order to allow for the 
magnetic calculation (either colinear or non-colinear).

## Projections
These projections will mainly be used for generating Wannier90 inputs, and to distinguish which indices
in the Wannier90 output matrices correspond to the various atoms/orbitals.
```@docs
Projection
```

### Orbitals
```@docs
Structures.Orbital
```

## DFT + U 
```@docs
DFTU
```
