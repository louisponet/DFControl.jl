# DFControl
[![Build Status](https://travis-ci.com/louisponet/DFControl.jl.svg?branch=master)](https://travis-ci.com/louisponet/DFControl.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/05vktbyj39u47usa?svg=true)](https://ci.appveyor.com/project/louisponet/dfcontrol-jl)
[![Coverage Status](https://coveralls.io/repos/github/louisponet/DFControl.jl/badge.svg?branch=master)](https://coveralls.io/github/louisponet/DFControl.jl?branch=master)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://louisponet.github.io/DFControl.jl/stable)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://louisponet.github.io/DFControl.jl/dev)

This package is a tool to interact with DFT related packages. Currently best support is for Quantum-Espresso and WANNIER90.
The support for Abinit is highly experimental and will get updated very soon.
There is some integration with Juno, namely the display of various Types is specifically tuned for ease of use.

## Installation

This package is registered, so as any normal Julia package:
```julia
Pkg.add("DFControl")
```

See the [documentation](https://louisponet.github.io/DFControl.jl/) on how to get started.
