const celldm_1 = Symbol("celldm(1)")

function read_block(f, startstr::String, endstr::String)
    block = [startstr]
    line = readline(f)
    while !eof(f)
        while !occursin(endstr, line)
            line = readline(f)
            push!(block, line)
        end
        return block
    end
    error("Block not found: start = $startstr, end = $endstr.")
end

#this is both flags and variables, QE calls it variables so ok
struct QEVariableInfo{T}
    name::Symbol
    typ::Type{T}
    # default::Union{T,Nothing,Symbol}  #check again v0.7 Some
    description::Array{String, 1}
end
QEVariableInfo{T}(name::Symbol, description) where T = QEVariableInfo{T}(name, T, description)
QEVariableInfo() = QEVariableInfo(:error, Nothing, String[])

function read_qe_variable(lines, i)
    name = gensym()
    var_i = i
    i += 2
    typ = fort2julia(strip_split(lines[i])[2])
    description = String[]
    # default = nothing
    i += 1
    line = lines[i]
    while !occursin("+------", line)
        # if occursin("Default", line)
        #     _t = strip_split(line)[2]
        #     _t = strip(strip(_t,'('),')')
        #     if occursin("D", _t)
        #         default = Meta.parse(typ, replace(_t,"D","e"))
        #     else
        #         _t = occursin("=",_t) ?split(_t,"=")[end] : _t
        #         default = typ ==String ? _t : Meta.parse(_t)
        #         println(_t)
        #         if typeof(default) != Symbol
        #             default = convert(typ, default)
        #         end
        #     end
        if occursin("Description", line)
            push!(description, strip_split(line,":")[2])
            i += 1
            line = lines[i]
            while !occursin("+------", line)
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
    line = lines[var_i]
    if occursin("Variables", line)
        spl = [split(x,"(")[1] for x in strip.(filter(x -> !occursin("=", x), split(line)[2:end]), ',')]
        names = Symbol.(spl)
        return [QEVariableInfo{typ}(name, description) for name in names], i
    else
        if occursin("(", line) && occursin(")", line)
            name = Symbol(split(strip_split(line, ":")[end],"(")[1])
        else
            name = Symbol(strip_split(line,":")[end])
        end
        return [QEVariableInfo{typ}(name, description)], i
    end
end

abstract type AbstractBlockInfo end

struct QEControlBlockInfo <: AbstractBlockInfo
    name::Symbol
    variables::Array{<:QEVariableInfo, 1}
end
function QEControlBlockInfo(lines::Array{<:AbstractString, 1})
    name  = Symbol(lowercase(strip_split(lines[1], "&")[2]))
    varinfos = QEVariableInfo[]
    for i=1:length(lines)
        line = lines[i]
        if occursin("Variable", line)
            variables, _i = read_qe_variable(lines, i)
            i += _i
            for var in variables
                push!(varinfos, var)
            end
        end
    end
    return QEControlBlockInfo(name, varinfos)
end

#TODO rewrite this this is garbage
function qevariable(block::AbstractBlockInfo, variable_name::Symbol)
    varstr1 = string(variable_name)
    for var in block.variables
        varstr2 = string(var.name)
        if occursin(varstr2, varstr1) && (length(varstr1) == length(varstr2) || length(varstr1) == length(varstr2) + 2)
            return var
        end
    end
    return QEVariableInfo()
end

struct QEDataBlockInfo <: AbstractBlockInfo
    name                ::Symbol
    description         ::Array{String, 1}
    options             ::Array{Symbol, 1}
    options_description ::Array{String, 1}
    variables           ::Array{<:QEVariableInfo, 1}
end
function QEDataBlockInfo(lines::Array{<:AbstractString, 1})
    spl                 = split(lines[1])
    name                = lowercase(spl[2])
    options             = Symbol.(spl[4:2:end])
    description         = String[]
    variables           = QEVariableInfo[]
    options_description = String[]
    i = 2
    while i <= length(lines) - 1
        line = strip(lines[i])
        if occursin("___________", line) && !occursin("+--", line)
            i += 1
            line = lines[i]
            while !occursin("----------------", line)
                push!(description, strip(line))
                i += 1
                line = lines[i]
            end
            i += 1

        elseif occursin("Card's flags:", line)
            i += 2
            if !occursin("Description:", lines[i])
                i += 1
            end
            push!(options_description,  join(split(lines[i])[2:end]," "))
            i += 1
            line = lines[i]
            while !occursin("----------------", line)
                push!(options_description, strip(line))
                i += 1
                line = lines[i]
            end
            i += 1
        elseif occursin("Variable", line)
            vars, _i = read_qe_variable(lines, i)
            i = _i
            for var in vars
                push!(variables, var)
            end
        end
        i += 1
    end
    return QEDataBlockInfo(Symbol(name), description, options, options_description, variables)
end

struct QEInputInfo
    exec::String
    control::Vector{QEControlBlockInfo}
    data::Vector{QEDataBlockInfo}
end

function QEInputInfo(filename::String; exec_name = join([lowercase(splitext(split(filename, "_")[end])[1]),".x"],""))
    control_block_infos = QEControlBlockInfo[]
    data_block_infos = QEDataBlockInfo[]
    open(filename, "r") do f
        while !eof(f)
            line = readline(f)
            if occursin("NAMELIST", line)
                push!(control_block_infos,
                      QEControlBlockInfo(read_block(f, line, "END OF NAMELIST")))
            elseif occursin("CARD:", line)
                push!(data_block_infos,
                      QEDataBlockInfo(read_block(f, line, "END OF CARD")))
            end
        end
    end
    return QEInputInfo(exec_name, control_block_infos, data_block_infos)
end

allflags(info::QEInputInfo) = flatten([[i.variables for i in info.control]; [i.variables for i in info.data]])

const QEInputInfos = begin
    input_files = searchdir(joinpath(@__DIR__, "../../assets/inputs/qe/"), "INPUT")
    file_paths  = joinpath(@__DIR__, "../../assets/inputs/qe/") .* input_files
    pw2wannier90_flags = [QEVariableInfo{String}(:outdir        , ["location of temporary output files"]),
                          QEVariableInfo{String}(:prefix        , ["pwscf filename prefix"]),
                          QEVariableInfo{String}(:seedname      , ["wannier90 input/output filename prefix"]),
                          QEVariableInfo{String}(:wan_mode      , ["'standalone' or 'library'"]),
                          QEVariableInfo{String}(:spin_component, ["'none', 'up' or 'down'"]),
                          QEVariableInfo{Bool}(:write_mmn     , ["compute M_mn matrix"]),
                          QEVariableInfo{Bool}(:write_amn     , ["compute A_mn matrix"]),
                          QEVariableInfo{Bool}(:write_unk     , ["write wavefunctions to file"]),
                          QEVariableInfo{Bool}(:write_uHu     , ["write the hamiltonian elements between different k-values"]),
                          QEVariableInfo{Bool}(:wvfn_formatted, ["formatted or unformatted output for wavefunctions"]),
                          QEVariableInfo{Bool}(:reduce_unk    , ["output wavefunctions on a coarse grid to save memory"])]
    pw2wannier90_block_info = QEControlBlockInfo(:inputpp, pw2wannier90_flags)
    vcat(QEInputInfo.(file_paths), QEInputInfo("pw2wannier90.x", [pw2wannier90_block_info], QEDataBlockInfo[]))
end

qe_input_info(input::DFInput{QE}) = getfirst(x-> occursin(x.exec, input.exec), QEInputInfos)
qe_input_info(exec::AbstractString) = getfirst(x-> occursin(x.exec, exec), QEInputInfos)
qe_input_flags(exec::AbstractString) = allflags(qe_input_info(exec))

function qevariable(input_info::QEInputInfo, variable_name::Symbol)
    for block in vcat(input_info.control, input_info.data)
        var = qevariable(block, variable_name)
        if var.typ != Nothing
            return var
        end
    end
    return QEVariableInfo()
end

function qevariable(variable_name::Symbol)
    for info in QEInputInfos
        var = qevariable(info, variable_name)
        if var.typ != Nothing
            return var
        end
    end
    return QEVariableInfo()
end

function qe_block_variable(input_info::QEInputInfo, variable_name)
    for block in vcat(input_info.control, input_info.data)
        var = qevariable(block, variable_name)
        if var.typ != Nothing
            return block, var
        end
    end
    return :error, QEVariableInfo()
end

function qevariable(exec::Exec, varname)
    for input_info in QEInputInfos
        if occursin(input_info.exec, exec.exec)
            return qevariable(input_info, varname)
        end
    end
    return QEVariableInfo()
end


function qe_block_info(block_name::Symbol)
    for input_info in QEInputInfos
        for block in [input_info.control;input_info.data]
            if block.name == block_name
                return block
            end
        end
    end
end


all_qe_block_flags(input::DFInput{QE}, block_name) = getfirst(x -> x.name == block, qe_input_info(input).control).variables
all_qe_block_flags(exec::AbstractString, block_name) = getfirst(x -> x.name == block_name, qe_input_info(exec).control).variables

function qe_block_variable(exec::AbstractString, flagname)
    for input_info in QEInputInfos
        if occursin(input_info.exec, exec)
            return qe_block_variable(input_info, flagname)
        end
    end
    return :error, QEVariableInfo()
end

qe_block_variable(input::DFInput, flagname) = qe_block_variable(execs(input)[2].exec, flagname)

flagtype(input::DFInput{QE}, flag) = qevariable(execs(input)[2], flag).typ
flagtype(::Type{QE}, exec, flag) = qevariable(exec, flag).typ
