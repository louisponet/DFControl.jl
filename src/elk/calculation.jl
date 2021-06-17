infilename(calculation::DFCalculation{Elk}) = "elk.in"

isbandscalc(calculation::DFCalculation{Elk})   = calculation.name == "20"

isnscfcalc(calculation::DFCalculation{Elk})    = calculation.name == "elk2wannier" #nscf == elk2wan??

isscfcalc(calculation::DFCalculation{Elk})     = calculation.name ∈ ["0", "1"]
