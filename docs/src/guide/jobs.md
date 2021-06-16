# Jobs

## Basics
A `DFJob` is a collection of `DFInputs` that can be scheduled to run or not in a linear fashion in the job's `local_dir`, using the `Structure` associated with the job. 

There are three ways of constructing a job:
- `DFJob(name::String, structure::Structure, inputs::Vector{DFInput}, input_flags::Pair{Symbol, Any}...; local_dir::String = "")` - creates a new job 
- `DFJob("/abs/path/to/job/dir")` - loads a job from the specified path
- `DFJob("<job dir fuzzy>")` - looks through the jobs known to `DFControl` (i.e. that were saved previously) and shows a menu allowing to select the desired one

A job can be saved through `save(job)`, this will generate and save all the input files associated with each `DFInput` in the job, and subsequently save a `job.tt` submission script. `job.local_dir` can be inspected to get an idea of what these are.

`submit` will first `save` the job and then attempt to run the job script, first using `sbatch job.tt` followed by `bash job.tt`. 

The job script preamble is defined by a `Vector{String}` in `job.header`, which will be pasted as lines before the execution lines. This can be useful for specifying `#SBATCH` options, `module load` directives, `export OMP_NUM_THREADS=1` environment variables, and similar.

To change a string inside one of the `header` lines, `set_headerword!(job, "string123" => "string124")` can be used, where any string that matches `"string123"` will be replaced by `"string124"`.


## 
