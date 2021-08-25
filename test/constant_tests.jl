using DFControl, Test

@test eltype(Calculations.qe_flaginfo(Exec(; exec = "projwfc.x"), :calculation)) == Nothing
@test eltype(Calculations.qe_flaginfo(Exec(; exec = "pw.x"), :calculation)) ==
      eltype(Calculations.qe_flaginfo(:calculation))
