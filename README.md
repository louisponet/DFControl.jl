# DFControl
[![Build Status](https://travis-ci.org/louisponet/DFControl.jl.svg?branch=master)](https://travis-ci.org/louisponet/DFControl.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/05vktbyj39u47usa?svg=true)](https://ci.appveyor.com/project/louisponet/dfcontrol-jl)
[![Coverage Status](https://coveralls.io/repos/github/louisponet/DFControl.jl/badge.svg?branch=master)](https://coveralls.io/github/louisponet/DFControl.jl?branch=master)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://louisponet.github.io/DFControl.jl/latest)

This package is a tool to interact with DFT related packages. Currently best support is for Quantum-Espresso and WANNIER90.
The support for Abinit is highly experimental and will get updated very soon.
There is some integration with Juno, namely the display of various Types is specifically tuned for ease of use.

## Installation

This package is registered, but the recommended way to install it is by:
```julia
Pkg.clone("https://github.com/louisponet/DFControl.jl.git")
Pkg.build("DFControl")
```
Since the package has changed quite a lot since the last release.

This will create a directory `user_defaults` with file `user_defaults.jl` inside the `DFControl` source folder. This is done because it allows one to define certain variables and defaults that will get loaded when `using DFControl` is called. The main use for this is to define various defaults, which make a lot of actions more streamlined.

## Defaults/Setup

After installation it is recommended to define a couple of those defaults.
```julia
using DFControl
```
default server:
```julia
setdefault_server("blabla@server.com") #default server
```

default pseudo potentials:
this defines a directory on the server to look through for the pseudo potentials for each element. This will be used for certain options when changing atom properties etc.
```julia
setdefault_pseudodir(:pbesol,  "pseudos/pbesol/") #change to your pseudo_set_name and directory of choice
setdefault_pseudodir(:pbesolrel, "pseudos/pbesolrel/")
#more sets can be defined
```

followed by the actual loading of the various pseudo filenames:
```julia
configuredefault_pseudos()
```
This will then connect to the server, look through all the defined pseudo sets inside the `default_pseudos` and tries to link for each element, for each set the correct filename.
To find out the filename of a certain atom for a certain pseudo set, or to check whether your config worked, you can do:
```julia
getdefault_pseudo(:O, :pbesolrel) #again change `pbesolrel` to the set you defined before
```
This should return you the filename for the pseudo potential file of Oxygen, `:O`, the format for elements in general is e.g. `:Mn`.
If multiple pseudos are defined for one set and element, you can specify keyword `pseudo_fuzzy = ...` to pull out the one you want to use.
For more info on other default functionality please look in the documentation and examples.

## General Usage

The main types in around which the package revolves are the `DFJob`, `DFInput`.
A `DFJob` is comprised of a `Structure`, representing the structure that is simulated in the job, and a collection of calculations to be done. Other fields in `DFJob` are auxiliary properties, such as name, directories, etc.

At this point, the main way the package works is by reading slurm job scripts such as the `test/test_job/job.tt` one.
What is most important for succesfully parsing these is the format used in the lines that do the actual calculations, i.e. `runcommand exec <input_file.in> output_file.out`, for example: `mpirun -np 24 ~/bin/projwfc.x  <projwfc.in> projwfc.out`.
The parsing is sort of robust but one should probably stick to this format.
All other not recognized lines will be saved in the `header` field of the `DFJob`.
The calculations in commented out lines will also be read and loaded, but they will have a field `DFInput.run=false`.
When the `job.tt` file gets written upon saving of the `DFJob`, calculations which are marked to not run will be written in the `job.tt` file, but be commented out.

As a quick start to see this in action you can do (Juno is highly recommended, for reading clarity)
```julia
job = DFJob(joinpath(dirname(pathof("DFControl")), "..", "test/test_job/"))
```

It will automatically look through the directory for a file which matches the fuzzy `*job*`. This can be specified through kwargs, further info in the documentation and examples.

To do something similar on a directory on a server, you can do
```julia
job = DFJob("path/to/job/starting_from_home", "local_dir", server=getdefault_server())
```
If the `local_dir` doesn't exist it will be first created before pulling the `job` script and the calculations that it read from this file.

If you ran some calculations, you should be able to pull the outputs into the `local_dir` of the `job`, followed by reading them:
```julia
outs = pulloutputs(job)
outdata= outputdata(job)
```

If one of the calculations that were performed was a calculation that produced `bands` of some sort (currently both outputs of `nscf` and `bands` calculations count) using Quantum-Espresso, you can do
```julia
bands = outdata["name_of_calculation(e.g. nscf)"]
using Plots
plot(bands, fermi=outdata["nscf"][:fermi])
```

Calculations can be set to run or not by
```julia
setflow!(job, "" => false) #none will run
setflow!(job, "nscf" => true, "bands" => true) #args..., and they are matched fuzzily
```

A job can be submitted by
```julia
submit(job)
```
If the job had "localhost" as it's server, it will run `qsub job.tt` locally, whereas if server and server_dir are something else, it will push the files and subsequently run them on the server.

This gave just a very small overview of the functionality, please look into the [documentation](https://louisponet.github.io/DFControl.jl/latest) and the examples for more.
