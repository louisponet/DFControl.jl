using DFControl, Test

@test Calculations.flagtype(Calculation{QE}(name = "dfd", exec = Exec(; path = "projwfc.x")), :calculation) == Nothing
@test Calculations.flagtype(Calculation{QE}(name = "dfd", exec = Exec(; path = "pw.x")), :calculation) == String
