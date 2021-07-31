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

ψ_cutoff_flag(c::DFCalculation) = ψ_cutoff_flag(package(c))
ρ_cutoff_flag(c::DFCalculation) = ρ_cutoff_flag(package(c))

function pdos(calculation::DFCalculation, args...)
    @error "pdos reading not implemented for package $(package(calculation))."
end

function Emin_from_projwfc(calculation::DFCalculation, args...)
    @error "Emin_from_projwfc is not implemented for package $(package(calculation))."
end

function readoutput(c::DFCalculation, args...; kwargs...)
    @error "Output parsing for package $(package(c)) not implemented."
end

rm_outfiles(calc::DFCalculation) = rm.(outfiles(calc))

include("qe/calculation.jl")
include("elk/calculation.jl")
include("wannier90/calculation.jl")
