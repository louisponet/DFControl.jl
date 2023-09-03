foundflags = Calculations.documentation(QE, "noncollinear")
@test length(foundflags) == 2
@test eltype(Calculations.documentation(QE, :electron_maxstep)) == Int
