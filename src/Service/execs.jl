using ..DFControl: RUN_EXECS, QE_EXECS, WAN_EXECS, ELK_EXECS

function calculationparser(exec::Exec)
    if DFC.is_qe_exec(exec)
        qe_read_calculation
    elseif DFC.is_wannier_exec(exec)
        wan_read_calculation
    elseif DFC.is_elk_exec(exec)
        elk_read_calculation
    end
end

allexecs() = vcat(RUN_EXECS, QE_EXECS, WAN_EXECS, ELK_EXECS)
parseable_execs() = vcat(QE_EXECS, WAN_EXECS, ELK_EXECS)
has_parseable_exec(l::String) = occursin(">", l) && any(occursin.(parseable_execs(), (l,)))

isparseable(exec::Exec) = exec.exec ∈ parseable_execs()

function verify_exec(e::Exec)
    if !isempty(e.dir) && ispath(DFC.path(e))
        return true
    else
        return Sys.which(e.exec) !== nothing
    end
end 
