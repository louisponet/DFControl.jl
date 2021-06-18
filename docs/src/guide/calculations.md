# Calculations

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

## Generating new calculations
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

