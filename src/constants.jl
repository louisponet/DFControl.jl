const abi_conversions = Dict{Symbol,Any}(:ev => 1 / 27.2113845,
                                         :ha => 1.0, 
                                         :ry => 0.5, 
                                         :ang => 1.889716164632,
                                         :bohr => 1.0,
                                         :au => 1.0,
                                         :angstr => 1.889716164632,
                                         :angstrom => 1.889716164632,
                                         :bohrs => 1.0,
                                         :hartree => 1.0,
                                         :hartrees => 1.0,
                                         :k => 1/3.1577513e5,
                                         :rydberg => 0.5,
                                         :rydbergs => 0.5,
                                         :t =>1/2.35e5,
                                         :tesla => 1/2.35e5)
                                         
const conversions = Dict{Symbol,Float64}(:bohr2ang => 0.529177)

abi_conversion(s::String) = abi_conversions[Symbol(s)]
convert_2abi(value, s::String) = value * abi_conversions[Symbol(s)]
convert_2abi(value, s::Symbol) = value * abi_conversions[s]

const qe_input_files = search_dir(joinpath(@__DIR__, "../assets/inputs/qe/"), "INPUT")

function f2julia(f_type)
    f_type = lowercase(f_type)
    if f_type == "real"
        return Float32
    elseif f_type == "real(kind=dp)"
        return Float64
    elseif f_type == "complex(kind=dp)"
        return Complex{Float64}
    elseif contains(f_type, "character")
        return String
    elseif f_type == "string"
        return String
    elseif f_type == "integer"
        return Int
    elseif f_type == "logical"
        return Bool
    end
end

"Reads all possible quantum espresso flags"
function read_qe_flags(filename)
    control_blocks = QEControlBlock[]
    open(filename, "r") do f
        while !eof(f)
            line = readline(f)

            if contains(line, "NAMELIST")
                name  = Symbol(lowercase(strip_split(line, "&")[2]))
                flags = Dict{Symbol,Any}()
                line  = readline(f)
                while !contains(line, "END OF NAMELIST")
                    if contains(line, "Variable")
                        readline(f)
                        value = f2julia(strip_split(readline(f))[2])
                        t_line = strip_split(line)
                        if t_line[1] == "Variables:"
                            flgs = Symbol.(strip.(t_line[2:end], ','))
                            for fl in flgs
                                flags[fl] = value
                            end
                        else
                            t_line = t_line[2]
                            if contains(t_line, "(") && contains(t_line, ")")
                                flag = Symbol(split(strip_split(t_line, ",")[1], "(")[1])
                            else
                                flag = Symbol(t_line)
                            end
                            flags[flag] = value
                        end
                        
                    end
                    line = readline(f)
                end
                push!(control_blocks, QEControlBlock(name, flags))
            end
        end
    end
    return control_blocks
end

const QEControlBlocks = vcat([read_qe_flags(joinpath(@__DIR__, "../assets/inputs/qe/") * file) for file in qe_input_files]...)

begin 
    flags = filter(x -> x.name == :inputpp, QEControlBlocks)[1].flags
    flags[:outdir]         = String
    flags[:prefix]         = String
    flags[:seedname]       = String
    flags[:wan_mode]       = String
    flags[:write_mmn]      = Bool
    flags[:write_amn]      = Bool
    flags[:write_unk]      = Bool
    flags[:wvfn_formatted] = Bool
    flags[:reduce_unk]     = Bool
    flags[:spin_component] = String
end

get_qe_flags(block) = filter(x -> x.name == block, QEControlBlocks)[1].flags

function get_qe_flag_block_type(flag)
    for block in QEControlBlocks
        if haskey(block.flags, flag)
            return block.name, block.flags[flag]
        end
    end
    return :error, Void
end

function read_wan_control_flags(filename::String)
    out = Dict{Symbol,Type}()
    open(filename, "r") do f
        while !eof(f)
            line = readline(f)
            if line == "" || line[1] == '!'
                continue
            else
                s_line    = split(line)
                flag      = Symbol(split(s_line[end],"(")[1])
                fl_type   = f2julia(strip(s_line[1],','))
                out[flag] = fl_type
            end
        end
    end
    return out
end

const WannierControlFlags = read_wan_control_flags(joinpath(@__DIR__, "../assets/inputs/wannier/input_flags.txt"))

get_wan_flag_type(flag) = haskey(WannierControlFlags, flag) ? WannierControlFlags[flag] : Void

function construct_abi_flags()
    out      = Dict{Symbol,Type}()
    open(joinpath(assets_dir,"inputs/abinit/input_variables.txt"), "r") do f
        while !eof(f)
            spl = split(readline(f))
            out[parse(spl[1])] = eval(parse(spl[2]))
        end
    end
    return out
end

const AbinitFlags    = construct_abi_flags()

get_abi_flag_type(flag) = haskey(AbinitFlags, flag) ? AbinitFlags[flag] : Void

