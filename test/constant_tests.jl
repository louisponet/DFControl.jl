using DFControl, Test

@test eltype(Calculations.qe_flaginfo(Exec(; path = "projwfc.x"), :calculation)) == Nothing
@test eltype(Calculations.qe_flaginfo(Exec(; path = "pw.x"), :calculation)) ==
      eltype(Calculations.qe_flaginfo(:calculation))
