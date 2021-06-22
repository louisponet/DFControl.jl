const QE_EXECS = ["pw.x", "projwfc.x", "pp.x", "ld1.x", "ph.x", "pw2wannier90.x", "hp.x"]

const QE_EXECFLAGS = ExecFlag[ExecFlag(:nk, "kpoint-pools", Int,
                                       "groups k-point parallelization into nk processor pools",
                                       0, 1),
                              ExecFlag(:ntg, "task-groups", Int, "FFT task groups", 0, 1),
                              ExecFlag(:ndiag, "diag", Int,
                                       "Number of processes for linear algebra", 0, 1),
                              ExecFlag(:ni, "images", Int,
                                       "Number of processes used for the images", 0, 1)]

qe_execflag(flag::AbstractString) = getfirst(x -> x.name == flag, QE_EXECFLAGS)
qe_execflag(flag::Symbol) = getfirst(x -> x.symbol == flag, QE_EXECFLAGS)

function parse_qe_execflags(line::Vector{<:AbstractString})
    flags = ExecFlag[]
    i = 1
    while i <= length(line)
        s = strip(line[i], '-')
        push!(flags, ExecFlag(qe_execflag(Symbol(s)), parse(Int, line[i+1])))
        i += 2
    end
    return flags
end

is_qe_exec(exec::Exec) = exec.exec ∈ QE_EXECS
