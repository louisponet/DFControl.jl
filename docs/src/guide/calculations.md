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
DFCalculation
InputData
```

## Basic interaction
```@docs
set_name!
data
set_data!
set_data_option!
set_kpoints!
```

## Flags
A big part of working with DFT calculations is specifying the various calculation flags.
Remembering all the names of the flags, where they belong, and what types they are expected to be,
is quite complicated and can lead to easily made mistakes like typos.
DFControl tries to catch these as it knows which flags are allowed for which calculations.
It will report when a flag can not be found for a given [`DFCalculation`](@ref Calculations),
and it will also try to convert a flag value to the expected type.

```@docs
Base.getindex(::DFCalculation, ::Symbol)
Base.setindex!(::DFCalculation, ::Symbol, ::Any)
set_flags!
rm_flags!
```

## Execs
```@docs
Exec
exec
execs
set_execdir!
set_execflags!
rm_execflags!
```

## [Generating new calculations](@id calculation_generation)
```@docs
gencalc_vcrelax(::DFCalculation{QE}, ::NTuple{6, Int}, ::Any...)
gencalc_scf(::DFCalculation{QE}, ::NTuple{6, Int}, ::Any...)
gencalc_bands(::DFCalculation{QE}, ::Vector{<:NTuple{4}}, ::Any...)
gencalc_nscf(::DFCalculation{QE}, kpoints::NTuple{3, Int}, ::Any...)
gencalc_projwfc(::DFCalculation{QE}, ::Real, ::Real, ::Real)
gencalc_wan
```

## Output 
```@docs
outputdata
readfermi
readbands
```

