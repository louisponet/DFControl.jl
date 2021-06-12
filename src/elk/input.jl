infilename(input::DFInput{Elk}) = "elk.in"

isbandscalc(input::DFInput{Elk})   = input.name == "20"

isnscfcalc(input::DFInput{Elk})    = input.name == "elk2wannier" #nscf == elk2wan??

isscfcalc(input::DFInput{Elk})     = input.name ∈ ["0", "1"]
