#these are all the control data, they hold the flags that guide the calculation
using ..DFControl: package, flagtype, flags
outfiles(c::DFCalculation) = filter(ispath, [outpath(c)])

function pdos(calculation::DFCalculation, args...)
    @error "pdos reading not implemented for package $(package(calculation))."
end

function Emin_from_projwfc(calculation::DFCalculation, args...)
    @error "Emin_from_projwfc is not implemented for package $(package(calculation))."
end

function readoutput(c::DFCalculation; kwargs...)
    @error "Output parsing for package $(package(c)) not implemented."
end

rm_outfiles(calc::DFCalculation) = rm.(outfiles(calc))

"""
    outputdata(calculation::DFCalculation; extra_parse_funcs=[], print=true, overwrite=true)

If an output file exists for `calculation` this will parse it and return a `Dict` with the parsed data.
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

outputdata(job["scf"], extra_parse_funcs = ["number of atoms/cell" => qe_parse_nat])
```
"""
function outputdata(calculation::DFCalculation;
                    extra_parse_funcs::Vector{<:Pair{String}} = Pair{String}[],
                    print = true, overwrite = true)
    if DFC.hasoutput(calculation)
        if !overwrite && !isempty(DFC.outdata(calculation))
            return DFC.outdata(calculation)
        else
            t = readoutput(calculation; parse_funcs = extra_parse_funcs)
            calculation.outdata = t === nothing ?
                                  parse_file(outpath(calculation), extra_parse_funcs) : t
            return calculation.outdata
        end
    end
    print &&
        (@warn "No output data or output file found for calculation: $(name(calculation)).")
    return DFC.SymAnyDict()
end

include("qe/calculation.jl")
include("elk/calculation.jl")
include("wannier90/calculation.jl")
