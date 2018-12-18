foundflags = search_documentation(QE, "noncollinear")
@test length(foundflags) == 2

show(foundflags[1])

@test eltype(documentation(QE, :electron_maxstep)) == Int
