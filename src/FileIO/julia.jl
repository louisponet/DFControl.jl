function julia_parse_calculation(file)
    if file isa IO || !occursin("\n", file)
        contents = read(file, String)
    else
        contents = file
    end
    return (flags = Dict{Symbol, Any}(:source => contents), data = InputData[], structure=nothing)
end

function readoutput(c::Calculation{Julia}, file; kwargs...)
    return parse_file(file, Pair{String}[]; kwargs...)
end

function Base.write(f::IO, calculation::Calculation{Julia}, structure)
    if haskey(calculation, :source)
        write(f, calculation[:source])
    elseif haskey(calculation, :script)
        write(f, read(calculation[:script]))
    end
end

infile_outfile_str(c::Calculation{Julia}) = "$(c.infile) > $(c.outfile)"
