#these are all the control data, they hold the flags that guide the calculation
using ..DFControl: package, flagtype, flags
outfiles(c::DFCalculation) = filter(ispath, [outpath(c)])

"Runs through all the set flags and checks if they are allowed and set to the correct value"
function convert_flags!(calculation::DFCalculation)
    for (flag, value) in flags(calculation)
        flagtype_ = flagtype(calculation, flag)
        if flagtype_ == Nothing
            @warn "Flag $flag was not found in allowed flags for exec $(execs(calculation)[2]). Removing flag."
            rm_flags!(calculation, flag)
            continue
        end
        if !(isa(value, flagtype_) || eltype(value) <: flagtype_)
            try
                if isbitstype(eltype(value))
                    if length(value) > 1
                        calculation[flag] = convert(flagtype_, value)
                    else
                        calculation[flag] = convert(eltype(flagtype_), value)
                    end
                else
                    calculation[flag] = convert.(flagtype_, value)
                end
            catch
                error("Input $(name(calculation)): Could not convert :$flag of value $value to the correct type ($flagtype_), please set it to the correct type.")
            end
        end
    end
end

#TODO implement abinit and wannier90
"""
    sanitize_flags!(calculation::DFCalculation, str::DFC.AbstractStructure)

Cleans up flags, i.e. remove flags that are not allowed and convert all
flags to the correct types.
Tries to correct common errors for different calculation types.
"""
function sanitize_flags!(calculation::DFCalculation, str::DFC.AbstractStructure)
    return convert_flags!(calculation)
end

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
