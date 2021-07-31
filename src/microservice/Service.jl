module Service
    using ..DFControl
    save_job(j::DFJob) = DFControl.save(j)
    get_job(p::AbstractString) = DFJob(p)
end
