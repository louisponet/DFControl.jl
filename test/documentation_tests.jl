foundflags = Calculations.documentation(QE, "noncollinear")
@test length(foundflags) == 2

show(foundflags[1])

@test eltype(Calculations.documentation(QE, :electron_maxstep)) == Int

@test length(Calculations.documentation(Elk, "magnetism")) == 2
@test Calculations.documentation(Elk, "magnetism")[1][2][1] == Calculations.documentation(Elk, :cmagz)
@test Calculations.documentation(Elk, :tasks) ==
      Calculations.getfirst(x -> x.name == :tasks, Calculations.ELK_CONTROLBLOCKS)
@test Calculations.documentation(Elk, "version")[1][1] ==
      Calculations.getfirst(x -> x.name == :tasks, Calculations.ELK_CONTROLBLOCKS)
