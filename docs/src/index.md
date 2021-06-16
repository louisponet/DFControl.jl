# DFControl
## Introduction
The goal of DFControl is to alleviate some of the tedious day to day busy-work that material's scientists using DFT 
codes encounter. It aids in the creation, execution, monitoring, management, and post-processing of DFT jobs. 
The core framework is code-agnostic, however, many features are so far only implemented for [Quantum-Espresso](https://quantum-espresso.org) and [Wannier90](http://www.wannier.org/). Rudimentary support for [ABINIT](https://www.abinit.org) and [ELK](https://elk.sourceforge.io) is also present.

## Philosophy
DFControl is structured to mimic the often adopted strategy of a linear submission script that specifies which input files 
are ran with which codes producing which outputs. Everything revolves around the [`DFJob`](@ref Jobs) struct that holds a collection of [`DFInput`](@ref Inputs) which will be ran sequentially in the job's directory. A `DFJob` is therefore linked to a given directory.
Using the structure of the script that is created upon saving or submitting a job, a `DFJob` can be easily reloaded at a later stage to continue the investigation where left off.

DFControl is aimed to be user friendly first and foremost, with features being added as time progresses without changing this core value. On the other hand, a core requirement has always been that power users that know all the ins and outs of the codes that they run should be able to have all the control they want without incurring friction from the library.
DFControl therefore should never be a barrier to code functionality, just a tool to increase efficiency of its users.

## Highlighted features
- Creation and submission of a simple self-consistent-field calculation starting from a structure Cif file in less than 10 lines of code (see the [Basic Tutorial](@ref)]). 
- Automatic validation and conversion of input flags
- Tracking of jobs for ease of continuation at later times
- Ease of input generation
- Automatic plotting of available results using a single command
- Input flag sanity checks
- Rudimentary job versioning to never lose previous job's results even if running in the same directory
- Fully human readable and transparent job directory structure 
