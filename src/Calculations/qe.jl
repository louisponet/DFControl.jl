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

include(joinpath(DFC.DEPS_DIR, "qeflags.jl"))
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

function qe_calculation_info(calculation::Calculation{QE})
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

function qe_all_block_flags(calculation::Calculation{QE}, block_name)
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

function qe_exec(calculation::Calculation{QE})
    exec = getfirst(x -> x.exec ∈ QE_EXECS, execs(calculation))
    if exec === nothing
        error("Calculation $calculation does not have a valid QE executable, please set it first.")
    end
    return exec
end

function qe_block_variable(calculation::Calculation, flagname)
    return qe_block_variable(qe_exec(calculation).exec, flagname)
end

function flagtype(calculation::Calculation{QE}, flag)
    return eltype(qe_flaginfo(qe_exec(calculation), flag))
end
flagtype(::Type{QE}, exec, flag) = eltype(qe_flaginfo(exec, flag))

isbands(c::Calculation{QE})   = c[:calculation] == "bands"
isnscf(c::Calculation{QE})    = c[:calculation] == "nscf"
isscf(c::Calculation{QE})     = c[:calculation] == "scf"
isvcrelax(c::Calculation{QE}) = c[:calculation] == "vc-relax"
isrelax(c::Calculation{QE})   = c[:calculation] == "relax"

function ispw(c::Calculation{QE})
    return isbands(c) || isnscf(c) || isscf(c) || isvcrelax(c) || isrelax(c)
end

issoc(c::Calculation{QE})     = c[:lspinorb] == true

function ismagnetic(c::Calculation{QE})
    return (hasflag(c, :nspin) && c[:nspin] > 0.0) ||
           (hasflag(c, :total_magnetization) && c[:total_magnetization] != 0.0)
end

function outfiles(c::Calculation{QE})
    files = [outpath(c)]
    for (is, fuzzies) in zip(("projwfc.x", "hp.x"), (("pdos",), ("Hubbard_parameters",)))
        if any(x -> x.exec == is, c.execs)
            for f in fuzzies
                append!(files, searchdir(c, f))
            end
        end
    end
    return filter(ispath, unique(files))
end

ψ_cutoff_flag(::Calculation{QE}) = :ecutwfc
ρ_cutoff_flag(::Calculation{QE}) = :ecutrho

for f in (:_cp, :_mv)
    base_func = Symbol(string(f)[2:end])
    @eval function $f(i::Calculation{QE}, dest::String; kwargs...)
        $f(inpath(i), joinpath(dest, i.infile); kwargs...)
        if hasoutfile(i)
            Base.$base_func(outpath(i), joinpath(dest, i.outfile); kwargs...)
        end
        if any(x -> x.exec == "projwfc.x", c.execs)
            for f in searchdir(i, "pdos")
                Base.$base_func(f, joinpath(dest, splitdir(f)[end]); kwargs...)
            end
        elseif any(x -> x.exec == "hp.x", c.execs)
            for f in searchdir(i, "Hubbard_parameters")
                Base.$base_func(f, joinpath(dest, splitdir(f)[end]); kwargs...)
            end
        end
        #TODO add ph.x outfiles
    end
end
