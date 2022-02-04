# Calculations

## Contents
```@contents
Pages=["calculations.md"]
```
## Index
```@index
Pages=["calculations.md"]
```

```@docs
Calculation
InputData
```

## Basic interaction
```@docs
Calculations.set_name!
Calculations.data
Calculations.set_kpoints!
```

## Flags
A big part of working with DFT calculations is specifying the various calculation flags.
Remembering all the names of the flags, where they belong, and what types they are expected to be,
is quite complicated and can lead to easily made mistakes like typos.
`DFControl` tries to catch these as it knows which flags are allowed for which calculations.
It will report when a flag can not be found for a given [`Calculation`](@ref Calculations),
and it will also try to convert a flag value to the expected type.

A [`Calculation`](@ref) behaves a lot as a normal `Dict` to interact with the stored flags.
This means that the usual `Dict` operations such as `haskey`, `get`, and `pop!` work on a
[`Calculation`](@ref).
```@docs
Base.getindex(::Calculation, ::Symbol)
Base.setindex!(::Calculation, ::Symbol, ::Any)
Calculations.set_flags!
```

## [Execs](@id execs_header)
```@docs
Exec
```

## [Generating new calculations](@id calculation_generation)
```@docs
Calculations.gencalc_vcrelax(::Calculation{QE}, ::NTuple{6, Int}, ::Any...)
Calculations.gencalc_scf(::Calculation{QE}, ::NTuple{6, Int}, ::Any...)
Calculations.gencalc_bands(::Calculation{QE}, ::Vector{<:NTuple{4}}, ::Any...)
Calculations.gencalc_nscf(::Calculation{QE}, kpoints::NTuple{3, Int}, ::Any...)
Calculations.gencalc_projwfc(::Calculation{QE}, ::Real, ::Real, ::Real)
Calculations.gencalc_wan
```

