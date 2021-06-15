# Installation

In case you don't have a working Julia installation yet, first
[download the Julia binaries](https://julialang.org/downloads/)
and follow the [Julia installation instructions](https://julialang.org/downloads/platform/).
DFControl is tested thoroughly with **Julia 1.6**, your mileage may vary with older versions.

Afterwards you can install DFControl
[like any other package](https://julialang.github.io/Pkg.jl/v1/getting-started/)
in Julia. For example run in your Julia REPL terminal:

```julia
import Pkg
Pkg.add("DFControl")
```
which will install the latest DFControl release.
Alternatively (if you like to be fully up to date) install the master branch:
```julia
import Pkg
Pkg.add(name="DFControl", rev="master")
```

DFControl is continuously tested on Debian, Ubuntu, mac OS and Windows and should work on
these operating systems out of the box.

After these steps it is highly recommended to go through some additional [Configuration](@ref).
