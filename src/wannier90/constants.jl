
const WannierControlFlags = Dict{Symbol, Type}()
function init_wan_control_flags()
    open(joinpath(@__DIR__, "..", "..", "assets", "inputs", "wannier", "input_flags.txt"), "r") do f
        while !eof(f)
            line = readline(f)
            if line == "" || line[1] == '!'
                continue
            else
                s_line    = split(line)
                flag      = Symbol(split(s_line[end],"(")[1])
                fl_type   = fort2julia(strip(s_line[1],','))
                WannierControlFlags[flag] = fl_type
            end
        end
    end
end
flagtype(::Type{Wannier90}, flag) = haskey(WannierControlFlags, flag) ? WannierControlFlags[flag] : Nothing
flagtype(::DFInput{Wannier90}, flag) = flagtype(Wannier90, flag)
