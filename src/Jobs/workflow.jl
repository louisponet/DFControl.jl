struct Workflow
    steps::Vector{Function}
    project_path::String
    modules::Vector{String}
end

function Workflow(steps::Function...)
    all = string.(values(Base.loaded_modules))
    valid = String[]
    ks = keys(Pkg.project().dependencies)
    for p in all
        if p ∈ ks
            push!(valid, p)
        end
    end
    return Workflow([steps...], Base.load_path()[1], valid)
end
