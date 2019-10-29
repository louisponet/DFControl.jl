const celldm_1 = Symbol("celldm(1)")
const QE_EXECS = [
    "pw.x",
    "projwfc.x",
    "pp.x",
    "ld1.x",
    "ph.x",
    "pw2wannier90.x"
]
#REVIEW: Should we make the flag name a String?
#QE calls these flags
struct QEFlagInfo{T}
    name::Symbol
    # default::Union{T,Nothing,Symbol}  #check again v0.7 Some
    description::String
end
QEFlagInfo() = QEFlagInfo{Nothing}(:error, "")
Base.eltype(x::QEFlagInfo{T}) where T = T

abstract type AbstractBlockInfo end

struct QEControlBlockInfo <: AbstractBlockInfo
    name::Symbol
    flags::Vector{<:QEFlagInfo}
end

#TODO rewrite this this is garbage
function qe_flaginfo(block::AbstractBlockInfo, variable_name::Symbol)
    varstr1 = string(variable_name)
    for var in block.flags
        varstr2 = string(var.name)
        if occursin(varstr2, varstr1) && (length(varstr1) == length(varstr2) || length(varstr1) == length(varstr2) + 2)
            return var
        end
    end
    return QEFlagInfo()
end

struct QEDataBlockInfo <: AbstractBlockInfo
    name                ::Symbol
    description         ::String
    options             ::Vector{Symbol}
    options_description ::String
    flags               ::Vector{<:QEFlagInfo}
end

struct QEInputInfo
    exec::String
    control::Vector{QEControlBlockInfo}
    data::Vector{QEDataBlockInfo}
end

allflags(info::QEInputInfo) = flatten([[i.flags for i in info.control]; [i.flags for i in info.data]])

include(joinpath(depsdir, "qeflags.jl"))
const QEInputInfos = _QEINPUTINFOS()
push!(QEInputInfos, QEInputInfo("pw2wannier90.x", [QEControlBlockInfo(:inputpp,[QEFlagInfo{String}(:outdir        , "location of temporary output files"),
                          QEFlagInfo{String}(:prefix        , "pwscf filename prefix"),
                          QEFlagInfo{String}(:seedname      , "wannier90 input/output filename prefix"),
                          QEFlagInfo{String}(:wan_mode      , "'standalone' or 'library'"),
                          QEFlagInfo{String}(:spin_component, "'none', 'up' or 'down'"),
                          QEFlagInfo{Bool}(:write_mmn     , "compute M_mn matrix"),
                          QEFlagInfo{Bool}(:write_amn     , "compute A_mn matrix"),
                          QEFlagInfo{Bool}(:write_unk     , "write wavefunctions to file"),
                          QEFlagInfo{Bool}(:write_uHu     , "write the hamiltonian elements between different k-values"),
                          QEFlagInfo{Bool}(:wvfn_formatted, "formatted or unformatted output for wavefunctions"),
                          QEFlagInfo{Bool}(:reduce_unk    , "output wavefunctions on a coarse grid to save memory")])], QEDataBlockInfo[]))

qe_input_info(input::DFInput{QE}) = getfirst(x-> occursin(x.exec, input.exec), QEInputInfos)
qe_input_info(exec::AbstractString) = getfirst(x-> occursin(x.exec, exec), QEInputInfos)
qe_input_flags(exec::AbstractString) = allflags(qe_input_info(exec))

function qe_flaginfo(input_info::QEInputInfo, variable_name::Symbol)
    for block in vcat(input_info.control, input_info.data)
        var = qe_flaginfo(block, variable_name)
        if eltype(var) != Nothing
            return var
        end
    end
    return QEFlagInfo()
end

function qe_flaginfo(variable_name::Symbol)
    for info in QEInputInfos
        var = qe_flaginfo(info, variable_name)
        if eltype(var) != Nothing
            return var
        end
    end
    return QEFlagInfo()
end

function qe_block_variable(input_info::QEInputInfo, variable_name)
    for block in vcat(input_info.control, input_info.data)
        var = qe_flaginfo(block, variable_name)
        if eltype(var) != Nothing
            return block, var
        end
    end
    return :error, QEFlagInfo()
end

function qe_flaginfo(exec::Exec, varname)
    for input_info in QEInputInfos
        if occursin(input_info.exec, exec.exec)
            return qe_flaginfo(input_info, varname)
        end
    end
    return QEFlagInfo()
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


qe_all_block_flags(input::DFInput{QE}, block_name) = getfirst(x -> x.name == block, qe_input_info(input).control).flags
qe_all_block_flags(exec::AbstractString, block_name) = getfirst(x -> x.name == block_name, qe_input_info(exec).control).flags

function qe_block_variable(exec::AbstractString, flagname)
    for input_info in QEInputInfos
        if occursin(input_info.exec, exec)
            return qe_block_variable(input_info, flagname)
        end
    end
    return :error, QEFlagInfo()
end

function qe_exec(input::DFInput{QE})
    exec = getfirst(x -> x.exec ∈ QE_EXECS, execs(input))
    if exec === nothing
        error("Input $input does not have a valid QE executable, please set it first.")
    end
    return exec
end

qe_block_variable(input::DFInput, flagname) = qe_block_variable(qe_exec(input).exec, flagname)

flagtype(input::DFInput{QE}, flag) = eltype(qe_flaginfo(qe_exec(input), flag))
flagtype(::Type{QE}, exec, flag) = eltype(qe_flaginfo(exec, flag))
