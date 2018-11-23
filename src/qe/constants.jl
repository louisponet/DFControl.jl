const celldm_1 = Symbol("celldm(1)")

#this is both flags and variables, QE calls it variables so ok
struct QEVariableInfo{T}
    name::Symbol
    typ::Type{T}
    # default::Union{T,Nothing,Symbol}  #check again v0.7 Some
    description::String
end
QEVariableInfo{T}(name::Symbol, description) where T = QEVariableInfo{T}(name, T, description)
QEVariableInfo() = QEVariableInfo(:error, Nothing, "")


abstract type AbstractBlockInfo end

struct QEControlBlockInfo <: AbstractBlockInfo
    name::Symbol
    variables::Vector{<:QEVariableInfo}
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
    description         ::String
    options             ::Vector{Symbol}
    options_description ::String
    variables           ::Vector{<:QEVariableInfo}
end

struct QEInputInfo
    exec::String
    control::Vector{QEControlBlockInfo}
    data::Vector{QEDataBlockInfo}
end

allflags(info::QEInputInfo) = flatten([[i.variables for i in info.control]; [i.variables for i in info.data]])

include(joinpath(depsdir, "qeflags.jl"))
const QEInputInfos = _QEINPUTINFOS()
push!(QEInputInfos, QEInputInfo("pw2wannier90.x", [QEControlBlockInfo(:inputpp,[QEVariableInfo{String}(:outdir        , "location of temporary output files"),
                          QEVariableInfo{String}(:prefix        , "pwscf filename prefix"),
                          QEVariableInfo{String}(:seedname      , "wannier90 input/output filename prefix"),
                          QEVariableInfo{String}(:wan_mode      , "'standalone' or 'library'"),
                          QEVariableInfo{String}(:spin_component, "'none', 'up' or 'down'"),
                          QEVariableInfo{Bool}(:write_mmn     , "compute M_mn matrix"),
                          QEVariableInfo{Bool}(:write_amn     , "compute A_mn matrix"),
                          QEVariableInfo{Bool}(:write_unk     , "write wavefunctions to file"),
                          QEVariableInfo{Bool}(:write_uHu     , "write the hamiltonian elements between different k-values"),
                          QEVariableInfo{Bool}(:wvfn_formatted, "formatted or unformatted output for wavefunctions"),
                          QEVariableInfo{Bool}(:reduce_unk    , "output wavefunctions on a coarse grid to save memory")])], QEDataBlockInfo[]))

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
