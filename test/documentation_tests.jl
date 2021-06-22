foundflags = documentation(QE, "noncollinear")
@test length(foundflags) == 2

show(foundflags[1])

@test eltype(documentation(QE, :electron_maxstep)) == Int

@test length(documentation(Elk, "magnetism")) == 2
@test documentation(Elk, "magnetism")[1][2][1] == documentation(Elk, :cmagz)
@test documentation(Elk, :tasks) ==
      DFControl.getfirst(x -> x.name == :tasks, DFControl.ELK_CONTROLBLOCKS)
@test documentation(Elk, "version")[1][1] ==
      DFControl.getfirst(x -> x.name == :tasks, DFControl.ELK_CONTROLBLOCKS)
