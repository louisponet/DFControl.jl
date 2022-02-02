#these are all the control data, they hold the flags that guide the calculation
function pdos(calculation::Calculation, args...)
    @error "pdos reading not implemented for package $(eltype(calculation))."
end

function Emin_from_projwfc(calculation::Calculation, args...)
    @error "Emin_from_projwfc is not implemented for package $(eltype(calculation))."
end

function readoutput(c::Calculation; kwargs...)
    @error "Output parsing for package $(eltype(c)) not implemented."
end

"""
    outputdata(calculation::Calculation, file; extra_parse_funcs=[], print=true, overwrite=true)

If `file` exists this will parse it and return a `Dict` with the parsed data.
If `overwrite=false` and `calculation.outputdata` is not empty, this will be returned instead of reparsing the
output file.

`extra_parse_funcs` should be `Vector{Pair{String,Function}}`, where the string will be used to match `occursin(str, line)`
for each line of the output file. If this returns `true`, the function will be called as `func(results_dict, line, file)`.
The purpose is to allow for additional parsing that is not implemented, or for temporary non-standard values that are printed
while working on the DFT code, e.g. for debugging.

!!! note
    This only works for files where not all information is pulled out for yet,
    e.g. projwfc.x outputs are fully parsed already.

Example (from src/qe/fileio.jl):
```julia

function qe_parse_nat(results, line, f)
    results[:nat] = parse(Int, split(line)[end])
end

outputdata(job["scf"],joinpath(job, job["scf"].outfile), extra_parse_funcs = ["number of atoms/cell" => qe_parse_nat])
```
"""
function outputdata(calculation::Calculation, file;
                    extra_parse_funcs::Vector{<:Pair{String}} = Pair{String}[],
                    print = true, overwrite = true)
    if ispath(file) 
        t = FileIO.readoutput(calculation, file; parse_funcs = extra_parse_funcs)
        return t === nothing ?
                              FileIO.parse_file(file,
                                                extra_parse_funcs) : t
    end
    print &&
        (@warn "File $file does not exist.")
    return Dict{Symbol,Any}()
end

verify_exec(args...) = Calculations.verify_exec(args...)

function save(e::Exec)
    @assert Calculations.isrunnable(e) "Exec is not runnable..."
    Calculations.save(e)
    return e
end
