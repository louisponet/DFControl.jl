cell(job::DFJob)         = cell(structure(job))
calculations(job::DFJob) = job.calculations
isarchived(job::DFJob) = occursin(".archived", job.local_dir)
