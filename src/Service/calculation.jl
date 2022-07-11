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

verify_exec(args...) = Calculations.verify(args...)
