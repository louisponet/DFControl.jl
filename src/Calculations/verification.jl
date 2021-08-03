const assets_dir = joinpath(@__DIR__, "..", "assets")
const conversions = Dict{Symbol,Float64}(:bohr2ang => 0.52917721092)
conversions[:ang2bohr] = 1 / conversions[:bohr2ang]
# include("abinit/constants.jl")
const celldm_1 = Symbol("celldm(1)")
#REVIEW: Should we make the flag name a String?

# QE
#QE calls these flags
struct QEFlagInfo{T}
    name::Symbol
    # default::Union{T,Nothing,Symbol}  #check again v0.7 Some
    description::String
end
QEFlagInfo() = QEFlagInfo{Nothing}(:error, "")
Base.eltype(x::QEFlagInfo{T}) where {T} = T

abstract type AbstractBlockInfo end

struct QEControlBlockInfo <: AbstractBlockInfo
    name::Symbol
    flags::Vector{<:QEFlagInfo}
end
Base.in(f::Symbol, i::QEControlBlockInfo) = findfirst(x -> x.name == f, i.flags) !== nothing


#TODO rewrite this this is garbage
function qe_flaginfo(block::AbstractBlockInfo, variable_name::Symbol)
    varstr1 = string(variable_name)
    for var in block.flags
        varstr2 = string(var.name)
        if occursin(varstr2, varstr1) &&
           (length(varstr1) == length(varstr2) || length(varstr1) == length(varstr2) + 2)
            return var
        end
    end
    return QEFlagInfo()
end

struct QEDataBlockInfo <: AbstractBlockInfo
    name                :: Symbol
    description         :: String
    options             :: Vector{Symbol}
    options_description :: String
    flags               :: Vector{<:QEFlagInfo}
end

struct QECalculationInfo
    exec::String
    control::Vector{QEControlBlockInfo}
    data::Vector{QEDataBlockInfo}
end

function allflags(info::QECalculationInfo)
    return flatten([[i.flags for i in info.control]; [i.flags for i in info.data]])
end

include(joinpath(depsdir, "qeflags.jl"))
const QECalculationInfos = _QEINPUTINFOS()
push!(QECalculationInfos,
      QECalculationInfo("pw2wannier90.x",
                        [QEControlBlockInfo(:inputpp,
                                            [QEFlagInfo{String}(:outdir,
                                                                "location of temporary output files"),
                                             QEFlagInfo{String}(:prefix,
                                                                "pwscf filename prefix"),
                                             QEFlagInfo{String}(:seedname,
                                                                "wannier90 calculation/output filename prefix"),
                                             QEFlagInfo{String}(:wan_mode,
                                                                "'standalone' or 'library'"),
                                             QEFlagInfo{String}(:spin_component,
                                                                "'none', 'up' or 'down'"),
                                             QEFlagInfo{Bool}(:write_spn,
                                                              "Write .spn matrix elements."),
                                             QEFlagInfo{Bool}(:write_mmn,
                                                              "compute M_mn matrix"),
                                             QEFlagInfo{Bool}(:write_amn,
                                                              "compute A_mn matrix"),
                                             QEFlagInfo{Bool}(:write_unk,
                                                              "write wavefunctions to file"),
                                             QEFlagInfo{Bool}(:write_uHu,
                                                              "write the hamiltonian elements between different k-values"),
                                             QEFlagInfo{Bool}(:wvfn_formatted,
                                                              "formatted or unformatted output for wavefunctions"),
                                             QEFlagInfo{Bool}(:reduce_unk,
                                                              "output wavefunctions on a coarse grid to save memory")])],
                        QEDataBlockInfo[]))

function qe_calculation_info(calculation::DFCalculation{QE})
    return getfirst(x -> occursin(x.exec, calculation.execs[end].exec), QECalculationInfos)
end
function qe_calculation_info(exec::AbstractString)
    return getfirst(x -> occursin(x.exec, exec), QECalculationInfos)
end
qe_calculation_flags(exec::AbstractString) = allflags(qe_calculation_info(exec))

function qe_flaginfo(calculation_info::QECalculationInfo, variable_name::Symbol)
    for block in vcat(calculation_info.control, calculation_info.data)
        var = qe_flaginfo(block, variable_name)
        if eltype(var) != Nothing
            return var
        end
    end
    return QEFlagInfo()
end

function qe_flaginfo(variable_name::Symbol)
    for info in QECalculationInfos
        var = qe_flaginfo(info, variable_name)
        if eltype(var) != Nothing
            return var
        end
    end
    return QEFlagInfo()
end

function qe_block_variable(calculation_info::QECalculationInfo, variable_name)
    for block in vcat(calculation_info.control, calculation_info.data)
        var = qe_flaginfo(block, variable_name)
        if eltype(var) != Nothing
            return block, var
        end
    end
    return :error, QEFlagInfo()
end

function qe_flaginfo(exec::Exec, varname)
    for calculation_info in QECalculationInfos
        if occursin(calculation_info.exec, exec.exec)
            return qe_flaginfo(calculation_info, varname)
        end
    end
    return QEFlagInfo()
end

function qe_block_info(block_name::Symbol)
    for calculation_info in QECalculationInfos
        for block in [calculation_info.control; calculation_info.data]
            if block.name == block_name
                return block
            end
        end
    end
end

function qe_all_block_flags(calculation::DFCalculation{QE}, block_name)
    return getfirst(x -> x.name == block, qe_calculation_info(calculation).control).flags
end
function qe_all_block_flags(exec::AbstractString, block_name)
    return getfirst(x -> x.name == block_name, qe_calculation_info(exec).control).flags
end

function qe_block_variable(exec::AbstractString, flagname)
    for calculation_info in QECalculationInfos
        if occursin(calculation_info.exec, exec)
            return qe_block_variable(calculation_info, flagname)
        end
    end
    return :error, QEFlagInfo()
end

function qe_exec(calculation::DFCalculation{QE})
    exec = getfirst(x -> x.exec ∈ QE_EXECS, execs(calculation))
    if exec === nothing
        error("Calculation $calculation does not have a valid QE executable, please set it first.")
    end
    return exec
end

function qe_block_variable(calculation::DFCalculation, flagname)
    return qe_block_variable(qe_exec(calculation).exec, flagname)
end

function flagtype(calculation::DFCalculation{QE}, flag)
    return eltype(qe_flaginfo(qe_exec(calculation), flag))
end
flagtype(::Type{QE}, exec, flag) = eltype(qe_flaginfo(exec, flag))

ψ_cutoff_flag(::DFCalculation{QE}) = :ecutwfc
ρ_cutoff_flag(::DFCalculation{QE}) = :ecutrho

# Wannier
include(joinpath(depsdir, "wannier90flags.jl"))
const WAN_FLAGS = _WAN_FLAGS()
flagtype(::Type{Wannier90}, flag) = haskey(WAN_FLAGS, flag) ? WAN_FLAGS[flag] : Nothing
flagtype(::DFCalculation{Wannier90}, flag) = flagtype(Wannier90, flag)


# ELK
include(joinpath(depsdir, "elkflags.jl"))

struct ElkFlagInfo{T}
    name::Symbol
    default::Union{T,Nothing}  #check again v0.7 Some
    description::String
end
ElkFlagInfo() = ElkFlagInfo{Nothing}(:error, "")
Base.eltype(x::ElkFlagInfo{T}) where {T} = T

struct ElkControlBlockInfo <: AbstractBlockInfo
    name::Symbol
    flags::Vector{<:ElkFlagInfo}
    description::String
end

const ELK_CONTROLBLOCKS = _ELK_CONTROLBLOCKS()

function elk_flaginfo(flag::Symbol)
    for b in ELK_CONTROLBLOCKS
        for f in b.flags
            if f.name == flag
                return f
            end
        end
    end
end

elk_block_info(name::Symbol) = getfirst(x -> x.name == name, ELK_CONTROLBLOCKS)

function elk_block_variable(flag_name::Symbol)
    for b in ELK_CONTROLBLOCKS
        for f in b.flags
            if f.name == flag_name
                return b
            end
        end
    end
end

flagtype(::DFCalculation{Elk}, flag::Symbol) = eltype(elk_flaginfo(flag))


# Abinit
# const ABI_UNIT_NAMES = [lowercase(s) for s in [
#         "au",
#         "Angstr", "Angstrom", "Angstroms", "Bohr", "Bohrs",
#         "eV", "Ha", "Hartree", "Hartrees", "K", "Ry", "Rydberg", "Rydbergs",
#         "T", "Tesla"]]
#
# const abi_conversions = Dict{Symbol,Any}(:ev => 1 / 27.2113845,
#                                          :ha => 1.0,
#                                          :ry => 0.5,
#                                          :ang => 1.889716164632,
#                                          :bohr => 1.0,
#                                          :au => 1.0,
#                                          :angstr => 1.889716164632,
#                                          :angstrom => 1.889716164632,
#                                          :bohrs => 1.0,
#                                          :hartree => 1.0,
#                                          :hartrees => 1.0,
#                                          :k => 1/3.1577513e5,
#                                          :rydberg => 0.5,
#                                          :rydbergs => 0.5,
#                                          :t =>1/2.35e5,
#                                          :tesla => 1/2.35e5)
#
# convert_2abi(value, s::String) = value * abi_conversions[Symbol(s)]
# convert_2abi(value, s::Symbol) = value * abi_conversions[s]
#
# function construct_abi_flags()
#     out      = Dict{Symbol, Type}()
#     open(joinpath(assets_dir,"calculations/abinit/calculation_variables.txt"), "r") do f
#         while !eof(f)
#             spl = split(readline(f))
#             out[Meta.parse(spl[1])] = eval(Meta.parse(spl[2]))
#         end
#     end
#     return out
# end
#

include(joinpath(depsdir, "abinitflags.jl"))
const AbinitFlags = _ABINITFLAGS()
#
# flagtype(calculation::DFCalculation{Abinit}, flag) = haskey(AbinitFlags, flag) ? AbinitFlags[flag] : Nothing
