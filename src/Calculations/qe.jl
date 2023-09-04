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

const QECalculationInfos = Ref(QECalculationInfo[])

function maybe_init_QECalculationInfos()
    if isempty(QECalculationInfos[])
        include(joinpath(DEPS_DIR, "qeflags.jl"))
        QECalculationInfos[] = Base.invokelatest(_QEINPUTINFOS,)
        push!(QECalculationInfos[],
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
    end
end

const QE7_2CalculationInfos = Ref(QECalculationInfo[])

function maybe_init_QE7_2CalculationInfos()
    if isempty(QE7_2CalculationInfos[])
        include(joinpath(DEPS_DIR, "qe7.2flags.jl"))
        QE7_2CalculationInfos[] = Base.invokelatest(_QE7_2INPUTINFOS, )
        push!(QE7_2CalculationInfos[],
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
    end
end

function qe_calculation_info(calculation::Calculation{QE})
    maybe_init_QECalculationInfos()        
    return getfirst(x -> occursin(x.exec, exec(calculation.exec)), QECalculationInfos[])
end

function qe_calculation_info(calculation::Calculation{QE7_2})
    maybe_init_QE7_2CalculationInfos()        
    return getfirst(x -> occursin(x.exec, exec(calculation.exec)), QE7_2CalculationInfos[])
end

function qe_calculation_info(exec::AbstractString)
    maybe_init_QECalculationInfos()        
    return getfirst(x -> occursin(x.exec, exec), QECalculationInfos[])
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
    maybe_init_QECalculationInfos()        
    for info in QECalculationInfos[]
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

function qe7_2_block_variable(exec::AbstractString, flagname)
    maybe_init_QE7_2CalculationInfos()        
    for calculation_info in QE7_2CalculationInfos[]
        if occursin(calculation_info.exec, exec)
            return qe_block_variable(calculation_info, flagname)
        end
    end
    return :error, QEFlagInfo()
end

function qe_block_variable(calculation::Calculation{QE}, flagname)
    e = exec(calculation.exec)
    if !(e ∈ QE_EXECS)
        error("Calculation $calculation does not have a valid QE executable, please set it first.")
    end
    
    maybe_init_QECalculationInfos()        
    for calculation_info in QECalculationInfos[]
        if occursin(calculation_info.exec, e)
            return qe_block_variable(calculation_info, flagname)
        end
    end
    return :error, QEFlagInfo()
    
end
function qe_block_variable(calculation::Calculation{QE7_2}, flagname)
    e = exec(calculation.exec)

    if !(e ∈ QE_EXECS)
        error("Calculation $calculation does not have a valid QE executable, please set it first.")
    end
    maybe_init_QE7_2CalculationInfos()        
    for calculation_info in QE7_2CalculationInfos[]
        if occursin(calculation_info.exec, e)
            return qe_block_variable(calculation_info, flagname)
        end
    end
    return :error, QEFlagInfo()
    
end

function flagtype(calculation::Calculation{QE}, flag)
    e = exec(calculation.exec)
    if !(e ∈ QE_EXECS)
        error("Calculation $calculation does not have a valid QE executable, please set it first.")
    end
    maybe_init_QECalculationInfos()        
    for calculation_info in QECalculationInfos[]
        if occursin(calculation_info.exec, e)
            return eltype(qe_flaginfo(calculation_info, flag))
        end
    end
    return eltype(QEFlagInfo())
end

function flagtype(calculation::Calculation{QE7_2}, flag)
    e = exec(calculation.exec)
    if !(e ∈ QE_EXECS)
        error("Calculation $calculation does not have a valid QE executable, please set it first.")
    end
    maybe_init_QE7_2CalculationInfos()        
    for calculation_info in QE7_2CalculationInfos[]
        if occursin(calculation_info.exec, e)
            return eltype(qe_flaginfo(calculation_info, flag))
        end
    end
    return eltype(QEFlagInfo())
end

isbands(c::Calculation{<:AbstractQE})   = get(c, :calculation, nothing) == "bands"
isnscf(c::Calculation{<:AbstractQE})    = get(c, :calculation, nothing) == "nscf"
isscf(c::Calculation{<:AbstractQE})     = get(c, :calculation, nothing) == "scf"
isvcrelax(c::Calculation{<:AbstractQE}) = get(c, :calculation, nothing) == "vc-relax"
isrelax(c::Calculation{<:AbstractQE})   = get(c, :calculation, nothing) == "relax"
isprojwfc(c::Calculation{<:AbstractQE}) = exec(c.exec) == "projwfc.x"
ishp(c::Calculation{<:AbstractQE})      = exec(c.exec) == "hp.x"

function ispw(c::Calculation{<:AbstractQE})
    return isbands(c) || isnscf(c) || isscf(c) || isvcrelax(c) || isrelax(c)
end

issoc(c::Calculation{<:AbstractQE}) = get(c, :lspinorb, false)

function ismagnetic(c::Calculation{<:AbstractQE})
    return get(c, :nspin, 0.0) > 0.0 || get(c, :total_magnetization, 0.0)
end

function outfiles(c::Calculation{<:AbstractQE})
    files = [c.outfile]
    for (is, fuzzies) in zip(("projwfc.x", "hp.x", "pp.x"), (("pdos",), ("Hubbard_parameters",), ("filplot", "fileout")))
        if c.exec.exec == is
            append!(files, fuzzies)
        end
    end
    return unique(files)
end

ψ_cutoff_flag(::Calculation{<:AbstractQE}) = :ecutwfc
ρ_cutoff_flag(::Calculation{<:AbstractQE}) = :ecutrho

function kgrid(na, nb, nc, ::Union{Type{QE}, Type{QE7_2}})
    return reshape([(a, b, c, 1 / (na * nb * nc))
                    for a in collect(range(0; stop = 1, length = na + 1))[1:end-1],
                        b in collect(range(0; stop = 1, length = nb + 1))[1:end-1],
                        c in collect(range(0; stop = 1, length = nc + 1))[1:end-1]],
                   (na * nb * nc))
end

function set_kpoints!(c::Calculation{<:AbstractQE}, k_grid::NTuple{3,Int}; print = true) #nscf
    print && !isnscf(c) && (@warn "Expected calculation to be 'nscf'.\nGot $c.")
    d = data(c, :k_points)
    if d !== nothing
        d.data = kgrid(k_grid..., c)
        d.option = :crystal
    else
        push!(c.data, InputData(:k_points, :crystal, kgrid(k_grid..., c)))
    end
    if prod(k_grid) > 100
        c[:control][:verbosity] = "high"
    end
    return c
end

function set_kpoints!(c::Calculation{<:AbstractQE}, k_grid::NTuple{6,Int}; print = true) #scf
    print &&
        !(isscf(c) || isvcrelax(c) || isrelax(c)) &&
        (@warn "Expected calculation to be scf, vc-relax, relax.")
    d = data(c, :k_points)
    if d !== nothing
        d.data = [k_grid...]
        d.option = :automatic
    else
        push!(c.data, InputData(:k_points, :automatic, [k_grid...]))
    end
    if prod(k_grid[1:3]) > 100
        c[:control][:verbosity] = "high"
    end
    return c
end

function set_kpoints!(c::Calculation{<:AbstractQE}, k_grid::Vector{<:NTuple{4}}; print = true,
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
        c[:control][:verbosity] = "high"
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
    gencalc_scf(template::Calculation{<:AbstractQE}, kpoints::NTuple{6, Int}, newflags...; name="scf")

Uses the information from the template and `supplied` kpoints to generate an scf calculation.
Extra flags can be supplied which will be set for the generated calculation.
"""
function gencalc_scf(template::Calculation{<:AbstractQE}, kpoints::NTuple{6,Int}, newflags...;
                     name = "scf")
    return calculation_from_kpoints(template, name, kpoints, :calculation => "scf",
                                    newflags...)
end

"""
    gencalc_vcrelax(template::Calculation{<:AbstractQE}, kpoints::NTuple{6, Int}, newflags...; name="scf")

Uses the information from the template and supplied `kpoints` to generate a vcrelax calculation.
Extra flags can be supplied which will be set for the generated calculation.
"""
function gencalc_vcrelax(template::Calculation{<:AbstractQE}, kpoints::NTuple{6,Int}, newflags...;
                         name = "vcrelax")
    return calculation_from_kpoints(template, name, kpoints, :calculation => "vc-relax",
                                    newflags...)
end

"""
    gencalc_bands(template::Calculation{<:AbstractQE}, kpoints::Vector{NTuple{4}}, newflags...; name="bands")

Uses the information from the template and supplied `kpoints` to generate a bands calculation.
Extra flags can be supplied which will be set for the generated calculation.
"""
function gencalc_bands(template::Calculation{<:AbstractQE}, kpoints::Vector{<:NTuple{4}}, newflags...;
                       name = "bands")
    return calculation_from_kpoints(template, name, kpoints, :calculation => "bands",
                                    newflags...)
end

"""
    gencalc_nscf(template::Calculation{<:AbstractQE}, kpoints::NTuple{3, Int}, newflags...; name="nscf")

Uses the information from the template and supplied `kpoints` to generate an nscf calculation.
Extra flags can be supplied which will be set for the generated calculation.
"""
function gencalc_nscf(template::Calculation{<:AbstractQE}, kpoints::NTuple{3,Int}, newflags...;
                      name = "nscf")
    return calculation_from_kpoints(template, name, kpoints, :calculation => "nscf",
                                    newflags...)
end

"""
    gencalc_projwfc(template::Calculation{<:AbstractQE}, Emin, Emax, DeltaE, newflags...; name="projwfc")

Uses the information from the template and supplied `kpoints` to generate a `projwfc.x` calculation.
Extra flags can be supplied which will be set for the generated calculation.
"""
function gencalc_projwfc(template::Calculation{<:AbstractQE}, Emin, Emax, DeltaE, extraflags...;
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
    exec = Exec(path=joinpath(dirname(template.exec), "projwfc.x"), modules = deepcopy(template.exec.modules))
    
    out = Calculation(deepcopy(template); name = name, exec = exec, data = InputData[])
    
    set_name!(out, "projwfc")
    
    empty!(out.flags)
    
    set_flags!(out, :Emin => Emin, :Emax => Emax, :DeltaE => DeltaE, :ngauss => ngauss,
               :degauss => degauss; print = false)
                
    set_flags!(out, extraflags...)
    return out
end
