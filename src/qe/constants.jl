#this is both flags and variables, QE calls it variables so ok
struct QEVariable{T}
    name::Symbol
    _type::Type{T}
    # default::Union{T,Void,Symbol}  #check again v0.7 Some
    description::Array{String, 1}
end
QEVariable{T}(name::Symbol, description) where T = QEVariable{T}(name, T, description)

function read_qe_variable(lines, i)
    name = gensym()
    var_i = i
    i += 2
    _type = fort2julia(strip_split(lines[i])[2])
    description = String[]
    # default = nothing
    i += 1
    line = lines[i]
    while !contains(line,"+------") 
        # if contains(line, "Default")
        #     _t = strip_split(line)[2]
        #     _t = strip(strip(_t,'('),')')
        #     if contains(_t, "D")
        #         default = parse(_type, replace(_t,"D","e"))
        #     else
        #         _t = contains(_t,"=") ?split(_t,"=")[end] : _t
        #         default = _type ==String ? _t : parse(_t)
        #         println(_t)
        #         if typeof(default) != Symbol 
        #             default = convert(_type, default)
        #         end
        #     end
        if contains(line, "Description")
            push!(description, strip_split(line,":")[2])
            i += 1
            line = lines[i]
            while !contains(line,"+------") 
                push!(description, strip(lines[i]))
                i += 1
                line = lines[i]
            end
            @goto break_label
        end
        i += 1
        line = lines[i]
    end
    @label break_label
    if contains(lines[var_i], "Variables")
        spl = [split(x,"(")[1] for x in strip.(filter(x -> !contains(x, "="), split(lines[var_i])[2:end]), ',')]
        names = Symbol.(spl)
        return [QEVariable{_type}(name, description) for name in names], i
    else
        if contains(lines[var_i], "(")&& contains(lines[var_i], ")")
            name = Symbol(split(strip_split(lines[var_i], ":")[end],"(")[1])
        else
            name = Symbol(strip_split(lines[var_i],":")[end])
        end
        return [QEVariable{_type}(name, description)], i
    end
end

struct QEControlBlockInfo
    name::Symbol
    flags::Array{<:QEVariable, 1}
end
function QEControlBlockInfo(lines::Array{<:AbstractString, 1})
    name  = Symbol(lowercase(strip_split(lines[1], "&")[2]))
    flags = QEVariable[] 
    for i=1:length(lines)
        line = lines[i]
        if contains(line, "Variable")
            variables, _i = read_qe_variable(lines, i)
            i += _i
            for var in variables
                push!(flags, var)
            end
        end
    end
    return QEControlBlockInfo(name, flags)
end

struct QEDataBlockInfo
    name                ::Symbol
    description         ::Array{String, 1}
    options             ::Array{Symbol, 1}
    options_description ::Array{String, 1}
    variables           ::Array{<:QEVariable, 1}
end
function QEDataBlockInfo(lines::Array{<:AbstractString, 1})
    spl                 = split(lines[1])
    name                = lowercase(spl[2])
    options             = Symbol.(spl[4:2:end])
    description         = String[]
    variables           = QEVariable[]
    options_description = String[]
    i = 2 
    while i <= length(lines) - 1
        line = strip(lines[i])
        if contains(line, "___________") && !contains(line, "+--")
            i += 1
            line = lines[i]
            while !contains(line, "----------------")
                push!(description, strip(line))
                i += 1
                line = lines[i]
            end
            i += 1

        elseif contains(line, "Card's flags:")
            i += 2
            if !contains(lines[i],"Description:")
                i += 1
            end
            push!(options_description,  join(split(lines[i])[2:end]," "))
            i += 1
            line = lines[i]
            while !contains(line, "----------------")
                push!(options_description, strip(line))
                i += 1
                line = lines[i]
            end
            i += 1
        elseif contains(line, "Variable")
            vars, _i = read_qe_variable(lines, i)
            i = _i
            for var in vars
                push!(variables, var)
            end
        end
        i += 1
    end
    return QEDataBlockInfo(name, description, options, options_description, variables)
end
struct QEInputInfo
    exec::String
    control_blocks::Array{QEControlBlockInfo, 1}
    data_blocks::Array{QEDataBlockInfo, 1}
end
function QEInputInfo(filename; exec_name = join([lowercase(splitext(split(filename, "_")[end])[1]),".x"],""))
    control_block_infos = QEControlBlockInfo[]
    data_block_infos = QEDataBlockInfo[]
    open(filename, "r") do f
        while !eof(f)
            line = readline(f)
            if contains(line, "NAMELIST")
                push!(control_block_infos,
                      QEControlBlockInfo(read_block(f, line, "END OF NAMELIST")))
            elseif contains(line, "CARD:")
                push!(data_block_infos,
                      QEDataBlockInfo(read_block(f, line, "END OF CARD")))
            end
        end
    end
    return QEInputInfo(exec_name, control_block_infos, data_block_infos)
end

const QEInputInfos = begin
    input_files = search_dir(joinpath(@__DIR__, "../../assets/inputs/qe/"), "INPUT")
    file_paths  = joinpath(@__DIR__, "../../assets/inputs/qe/") .* input_files
    pw2wannier90_flags = [QEVariable{String}(:outdir        , ["location of temporary output files"]),
                          QEVariable{String}(:prefix        , ["pwscf filename prefix"]),
                          QEVariable{String}(:seedname      , ["wannier90 input/output filename prefix"]),
                          QEVariable{String}(:wan_mode      , ["'standalone' or 'library'"]),
                          QEVariable{String}(:spin_component, ["'none', 'up' or 'down'"]),
                          QEVariable{Bool}(:write_mmn     , ["compute M_mn matrix"]),
                          QEVariable{Bool}(:write_amn     , ["compute A_mn matrix"]),
                          QEVariable{Bool}(:write_unk     , ["write wavefunctions to file"]),
                          QEVariable{Bool}(:wvfn_formatted, ["formatted or unformatted output for wavefunctions"]),
                          QEVariable{Bool}(:reduce_unk    , ["output wavefunctions on a coarse grid to save memory"])]
    pw2wannier90_block_info = QEControlBlockInfo(:inputpp, pw2wannier90_flags)
    vcat(QEInputInfo.(file_paths), QEInputInfo("pw2wannier90.x", [pw2wannier90_block_info], QEDataBlockInfo[]))
end

QEInputInfo(input::QEInput) = filter(x-> contains(input.exec, x.exec), QEInputInfos)[1]
QEInputInfo(exec::String) = filter(x-> contains(exec, x.exec), QEInputInfos)[1]

function get_qe_variable(input_info::QEInputInfo, variable_name)
    for block in vcat(input_info.control_blocks, input_info.data_blocks)
        for var in block.variables
            if var.name == variable_name
                return var
            end
        end
    end
end

function get_qe_block_variable(input_info::QEInputInfo, variable_name)
    for block in vcat(input_info.control_blocks, input_info.data_blocks)
        for var in block.variables
            if var.name == variable_name
                return block, var
            end
        end
    end
end

function get_qe_variable_info(input::QEInput, varname)
    for input_info in QEInputInfos
        if input_info.exec == input.exec
            return get_qe_variable(input_info, varname)
        end
    end
end

all_qe_block_flags(input::QEInput, block_name) = filter(x -> x.name == block, QEInputInfo(input).control_blocks)[1].flags
all_qe_block_flags(exec::String, block_name) = filter(x -> x.name == block_name, QEInputInfo(exec).control_blocks)[1].flags

function get_qe_flag_block_type(input::DFInput, flagname)
    for input_info in QEInputInfos
        if input_info.exec == input.exec
            return get_qe_block_variable(input_info, flagname)
        end
    end
    return :error, Void
end
