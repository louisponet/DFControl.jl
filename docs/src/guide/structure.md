# Structure
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
DFControl.Element
Atom
atoms(::Structure)
```
