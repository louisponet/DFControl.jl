function read_wan_control_flags(filename::String)
    out = OrderedDict{Symbol,Type}()
    open(filename, "r") do f
        while !eof(f)
            line = readline(f)
            if line == "" || line[1] == '!'
                continue
            else
                s_line    = split(line)
                flag      = Symbol(split(s_line[end],"(")[1])
                fl_type   = fort2julia(strip(s_line[1],','))
                out[flag] = fl_type
            end
        end
    end
    return out
end

const WannierControlFlags = read_wan_control_flags(joinpath(@__DIR__, "../../assets/inputs/wannier/input_flags.txt"))

get_wan_flag_type(flag) = haskey(WannierControlFlags, flag) ? WannierControlFlags[flag] : Void



