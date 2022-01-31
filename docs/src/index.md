# DFControl
## Introduction
The goal of DFControl is to alleviate some of the tedious day to day busy-work that material's scientists using DFT 
codes encounter. It aids in the creation, execution, monitoring, management, and post-processing of DFT jobs. 
The core framework is code-agnostic, however, most of the convenience features are so far only implemented for [Quantum-Espresso](https://quantum-espresso.org) and [Wannier90](http://www.wannier.org/). [ABINIT](https://www.abinit.org) and [ELK](https://elk.sourceforge.io) support is highly experimental and incomplete.

## Philosophy
DFControl is aimed to be user friendly first and foremost, with features being added as time progresses without changing this core value. On the other hand, a core requirement has always been that power users that know all the ins and outs of the codes that they run should be able to have all the control they want without incurring friction from the library.
DFControl therefore should never be a barrier to DFT code functionality, just a tool to increase efficiency of its users.

In light of this DFControl is structured to mimic the often adopted strategy of a linear submission script that specifies which input files 
are ran with which DFT codes producing which outputs. Everything revolves around the [`Job`](@ref Jobs) that holds a collection of [`Calculations`](@ref Calculations) which will be ran sequentially in the job's directory, on the job's crystal [`Structure`](@ref Structure). A [`Job`](@ref Jobs) is therefore identified with a specific directory.
Using the script that is created upon saving or submitting a job, a [`Job`](@ref Jobs) can be easily reloaded at a later stage to continue the investigation where left off.


## Highlighted features
- Creation and submission of a simple self-consistent-field calculation starting from a structure Cif file in less than 10 lines of code (see the [Basic Tutorial](@ref)])
- Separation of [`Structure`](@ref Structure) and [`Calculation`](@ref Calculations) so that the same [`Calculations`](@ref Calculations) can be copied and used with a different [`Structure`](@ref Structure) in a new job 
- Automatic validation and conversion of input flags
- [Tracking](@ref Registry) of [`Jobs`](@ref Job) for ease of continuation at later times
- Ease of [calculation generation](@ref calculation_generation)
- Automatic [plotting of available results](@ref results_plotting) using a single command
- Input flag sanity checks
- Rudimentary [`Job`](@ref) [versioning](@ref Versioning) to never lose previous [`Job's`](@ref Job) results even if running in the same directory
- Fully human readable and transparent [`Job`](@ref) directory structure 
