# DFControl
[![Build Status](https://travis-ci.org/louisponet/DFControl.jl.svg?branch=master)](https://travis-ci.org/louisponet/DFControl.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/05vktbyj39u47usa?svg=true)](https://ci.appveyor.com/project/louisponet/dfcontrol-jl)
[![Coverage Status](https://coveralls.io/repos/github/louisponet/DFControl.jl/badge.svg?branch=master)](https://coveralls.io/github/louisponet/DFControl.jl?branch=master)

This package is a tool to interact with DFT related packages. Currently best support is for Quantum-Espresso and WANNIER90, Abinit is also supported but limited.

## Installation

Since this package is not registered yet, the way to install it is by:
```julia
Pkg.clone("https://github.com/louisponet/DFControl.jl.git")
Pkg.build("DFControl")
```

This will create a directory `user_defaults` with file `user_defaults.jl` inside the `DFControl` source folder. This is done because it allows one to define certain variables and defaults that will get loaded when `using DFControl` is called. The main use for this is to define various defaults, which make a lot of actions more streamlined.

## Defaults/Setup

After installation it is recommended to define a couple of those defaults.

default server:
```julia
set_default_server("blabla@server.com") #default server
```

default pseudo potentials:
this defines a directory on the server to look through for the pseudo potentials for each element. This will be used for certain options when changing atom properties etc.
```julia
set_default_pseudo_dir(:pbesol,  "pseudos/pbesol/") #change to your pseudo_set_name and directory of choice
set_default_pseudo_dir(:pbesolrel, "pseudos/pbesolrel/")
#more sets can be defined
```

followed by the actual loading of the various pseudo filenames:
```julia
configure_default_pseudos()
```
This will then connect to the server, look through all the defined pseudo sets inside the `default_pseudos` and tries to link for each element, for each set the correct filename.
To find out the filename of a certain atom for a certain pseudo set, or to check whether your config worked, you can do:
```julia
get_default_pseudo(:O, :pbesolrel) #again change `pbesolrel` to the set you defined before
```
This should return you the filename for the pseudo potential file of Oxygen, `:O`, the format for elements in general is e.g. `:Mn`.
If multiple pseudos are defined for one set and element, you can specify keyword `pseudo_fuzzy = ...` to pull out the one you want to use.
For more info on other default functionality please look in the documentation and examples.

## General Usage
TODO
