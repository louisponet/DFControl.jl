const WAN_EXECS = ["wannier90.x"]

const WAN_EXECFLAGS = ExecFlag[ExecFlag(:pp, "preprocess", Nothing,
                                        "Whether or not to preprocess the wannier calculation",
                                        nothing, 1),]

wan_execflag(flag::AbstractString) = getfirst(x -> x.name == flag, WAN_EXECFLAGS)
wan_execflag(flag::Symbol) = getfirst(x -> x.symbol == flag, WAN_EXECFLAGS)
function parse_wan_execflags(line::Vector{<:AbstractString})
    flags = ExecFlag[]
    i = 1
    while i <= length(line)
        s = strip(line[i], '-')
        push!(flags, ExecFlag(wan_execflag(Symbol(s)), nothing))
        i += 1
    end
    return flags
end

is_wannier_exec(exec::Exec) = exec.exec ∈ WAN_EXECS
