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
    return vcat([[i.flags for i in info.control]; [i.flags for i in info.data]]...)
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
    return getfirst(x -> occursin(x.exec, calculation.exec.exec), QECalculationInfos)
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
    if !(calculation.exec.exec ∈ QE_EXECS)
        error("Calculation $calculation does not have a valid QE executable, please set it first.")
    end
    return calculation.exec
end

function qe_block_variable(calculation::Calculation, flagname)
    return qe_block_variable(qe_exec(calculation).exec, flagname)
end

function flagtype(calculation::Calculation{QE}, flag)
    return eltype(qe_flaginfo(qe_exec(calculation), flag))
end
flagtype(::Type{QE}, exec, flag) = eltype(qe_flaginfo(exec, flag))

isbands(c::Calculation{QE})   = get(c, :calculation, nothing) == "bands"
isnscf(c::Calculation{QE})    = get(c, :calculation, nothing) == "nscf"
isscf(c::Calculation{QE})     = get(c, :calculation, nothing) == "scf"
isvcrelax(c::Calculation{QE}) = get(c, :calculation, nothing) == "vc-relax"
isrelax(c::Calculation{QE})   = get(c, :calculation, nothing) == "relax"
isprojwfc(c::Calculation{QE}) = c.exec.exec == "projwfc.x"
ishp(c::Calculation{QE}) = c.exec.exec == "hp.x"

function ispw(c::Calculation{QE})
    return isbands(c) || isnscf(c) || isscf(c) || isvcrelax(c) || isrelax(c)
end

issoc(c::Calculation{QE}) = get(c, :lspinorb, false)

function ismagnetic(c::Calculation{QE})
    return get(c, :nspin, 0.0) > 0.0 || get(c, :total_magnetization, 0.0)
end

function outfiles(c::Calculation{QE})
    files = [c.outfile]
    for (is, fuzzies) in zip(("projwfc.x", "hp.x"), (("pdos",), ("Hubbard_parameters",)))
        if c.exec.exec == is
            append!(files, fuzzies)
        end
    end
    return unique(files)
end

ψ_cutoff_flag(::Calculation{QE}) = :ecutwfc
ρ_cutoff_flag(::Calculation{QE}) = :ecutrho

function kgrid(na, nb, nc, ::Type{QE})
    return reshape([(a, b, c, 1 / (na * nb * nc))
                    for a in collect(range(0; stop = 1, length = na + 1))[1:end-1],
                        b in collect(range(0; stop = 1, length = nb + 1))[1:end-1],
                        c in collect(range(0; stop = 1, length = nc + 1))[1:end-1]],
                   (na * nb * nc))
end

function set_kpoints!(c::Calculation{QE}, k_grid::NTuple{3,Int}; print = true) #nscf
    print && !isnscf(c) && (@warn "Expected calculation to be 'nscf'.\nGot $c.")
    d = data(c, :k_points)
    if d !== nothing
        d.data = kgrid(k_grid..., c)
        d.option = :crystal
    else
        push!(c.data, InputData(:k_points, :crystal, kgrid(k_grid..., c)))
    end
    prod(k_grid) > 100 && set_flags!(c, :verbosity => "high"; print = print)
    return c
end

function set_kpoints!(c::Calculation{QE}, k_grid::NTuple{6,Int}; print = true) #scf
    print &&
        !(isscf(c) || isvcrelax(c) || isrelax(c)) &&
        (@warn "Expected calculation to be scf, vc-relax, relax.\nGot $calc.")
    d = data(c, :k_points)
    if d !== nothing
        d.data = [k_grid...]
        d.option = :automatic
    else
        push!(c.data, InputData(:k_points, :automatic, [k_grid...]))
    end
    prod(k_grid[1:3]) > 100 && set_flags!(c, :verbosity => "high"; print = print)
    return c
end

function set_kpoints!(c::Calculation{QE}, k_grid::Vector{<:NTuple{4}}; print = true,
                      k_option = :crystal_b)
    print &&
        isbands(c) != "bands" &&
        (@warn "Expected calculation to be bands, got $(c[:calculation]).")
    @assert in(k_option, [:tpiba_b, :crystal_b, :tpiba_c, :crystal_c]) error("Only $([:tpiba_b, :crystal_b, :tpiba_c, :crystal_c]...) are allowed as a k_option, got $k_option.")
    if k_option in [:tpiba_c, :crystal_c]
        @assert length(k_grid) == 3 error("If $([:tpiba_c, :crystal_c]...) is selected the length of the k_points needs to be 3, got length: $(length(k_grid)).")
    end
    num_k = 0.0
    for k in k_grid
        num_k += k[4]
    end
    if num_k > 100.0
        set_flags!(c, :verbosity => "high"; print = print)
        if print
            @info "Verbosity is set to high because num_kpoints > 100,\n
                       otherwise bands won't get printed."
        end
    end
    d = data(c, :k_points)
    if d !== nothing
        data(c, :k_points).data = k_grid
        data(c, :k_points).option = k_option
    else
        push!(c.data, InputData(:k_points, k_option, k_grid))
    end
    return c
end

"""
    gencalc_scf(template::Calculation{QE}, kpoints::NTuple{6, Int}, newflags...; name="scf")

Uses the information from the template and `supplied` kpoints to generate an scf calculation.
Extra flags can be supplied which will be set for the generated calculation.
"""
function gencalc_scf(template::Calculation{QE}, kpoints::NTuple{6,Int}, newflags...;
                     name = "scf")
    return calculation_from_kpoints(template, name, kpoints, :calculation => "scf",
                                    newflags...)
end

"""
    gencalc_vcrelax(template::Calculation{QE}, kpoints::NTuple{6, Int}, newflags...; name="scf")

Uses the information from the template and supplied `kpoints` to generate a vcrelax calculation.
Extra flags can be supplied which will be set for the generated calculation.
"""
function gencalc_vcrelax(template::Calculation{QE}, kpoints::NTuple{6,Int}, newflags...;
                         name = "vcrelax")
    return calculation_from_kpoints(template, name, kpoints, :calculation => "vc-relax",
                                    newflags...)
end

"""
    gencalc_bands(template::Calculation{QE}, kpoints::Vector{NTuple{4}}, newflags...; name="bands")

Uses the information from the template and supplied `kpoints` to generate a bands calculation.
Extra flags can be supplied which will be set for the generated calculation.
"""
function gencalc_bands(template::Calculation{QE}, kpoints::Vector{<:NTuple{4}}, newflags...;
                       name = "bands")
    return calculation_from_kpoints(template, name, kpoints, :calculation => "bands",
                                    newflags...)
end

"""
    gencalc_nscf(template::Calculation{QE}, kpoints::NTuple{3, Int}, newflags...; name="nscf")

Uses the information from the template and supplied `kpoints` to generate an nscf calculation.
Extra flags can be supplied which will be set for the generated calculation.
"""
function gencalc_nscf(template::Calculation{QE}, kpoints::NTuple{3,Int}, newflags...;
                      name = "nscf")
    return calculation_from_kpoints(template, name, kpoints, :calculation => "nscf",
                                    newflags...)
end

"""
    gencalc_projwfc(template::Calculation{QE}, Emin, Emax, DeltaE, newflags...; name="projwfc")

Uses the information from the template and supplied `kpoints` to generate a `projwfc.x` calculation.
Extra flags can be supplied which will be set for the generated calculation.
"""
function gencalc_projwfc(template::Calculation{QE}, Emin, Emax, DeltaE, extraflags...;
                         name = "projwfc")
    occflag = get(template, :occupations, "fixed")
    ngauss  = 0
    if occflag == "smearing"
        smearingflag = get(template, :smearing,"gaussian")
        if smearingflag ∈ ("methfessel-paxton", "m-p", "mp")
            ngauss = 1
        elseif smearingflag ∈ ("marzari-vanderbilt", "cold", "m-v", "mv")
            ngauss = -1
        elseif smearingflag ∈ ("fermi-dirac", "f-d", "fd")
            ngauss = -99
        end
    end
    tdegaussflag = get(template, :degauss, nothing)
    degauss = tdegaussflag !== nothing ? tdegaussflag : 0.0
    exec = Exec(dir = template.exec.dir, exec = "projwfc.x", modules = template.exec.modules, flags = template.exec.flags)
    empty!(exec.flags)
    out = Calculation(deepcopy(template); name = name, exec = exec, data = InputData[])
    set_name!(out, "projwfc")
    empty!(out.flags)
    set_flags!(out, :Emin => Emin, :Emax => Emax, :DeltaE => DeltaE, :ngauss => ngauss,
               :degauss => degauss; print = false)
    set_flags!(out, extraflags...)
    return out
end
