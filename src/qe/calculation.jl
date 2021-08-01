function kgrid(na, nb, nc, ::Type{QE})
    return reshape([(a, b, c, 1 / (na * nb * nc))
                    for a in collect(range(0; stop = 1, length = na + 1))[1:end-1],
                        b in collect(range(0; stop = 1, length = nb + 1))[1:end-1],
                        c in collect(range(0; stop = 1, length = nc + 1))[1:end-1]],
                   (na * nb * nc))
end

function set_kpoints!(c::DFCalculation{QE}, k_grid::NTuple{3,Int}; print = true) #nscf
    calc = flag(c, :calculation)
    print && calc != "nscf" && (@warn "Expected calculation to be 'nscf'.\nGot $calc.")
    set_data!(c, :k_points, kgrid(k_grid..., c); option = :crystal, print = print)
    prod(k_grid) > 100 && set_flags!(c, :verbosity => "high"; print = print)
    return c
end

function set_kpoints!(c::DFCalculation{QE}, k_grid::NTuple{6,Int}; print = true) #scf
    calc = flag(c, :calculation)
    print &&
        calc != "scf" &&
        !occursin("relax", calc) &&
        (@warn "Expected calculation to be scf, vc-relax, relax.\nGot $calc.")
    set_data!(c, :k_points, [k_grid...]; option = :automatic, print = print)
    prod(k_grid[1:3]) > 100 && set_flags!(c, :verbosity => "high"; print = print)
    return c
end

function set_kpoints!(c::DFCalculation{QE}, k_grid::Vector{<:NTuple{4}}; print = true,
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
    set_data!(c, :k_points, k_grid; option = k_option, print = print)
    return c
end

"""
    gencalc_scf(template::DFCalculation{QE}, kpoints::NTuple{6, Int}, newflags...; name="scf")

Uses the information from the template and `supplied` kpoints to generate an scf calculation.
Extra flags can be supplied which will be set for the generated calculation.
"""
function gencalc_scf(template::DFCalculation{QE}, kpoints::NTuple{6,Int}, newflags...;
                     name = "scf")
    return calculation_from_kpoints(template, name, kpoints, :calculation => "scf",
                                    newflags...)
end

"""
    gencalc_vcrelax(template::DFCalculation{QE}, kpoints::NTuple{6, Int}, newflags...; name="scf")

Uses the information from the template and supplied `kpoints` to generate a vcrelax calculation.
Extra flags can be supplied which will be set for the generated calculation.
"""
function gencalc_vcrelax(template::DFCalculation{QE}, kpoints::NTuple{6,Int}, newflags...;
                         name = "vcrelax")
    return calculation_from_kpoints(template, name, kpoints, :calculation => "vc-relax",
                                    newflags...)
end

"""
    gencalc_bands(template::DFCalculation{QE}, kpoints::Vector{NTuple{4}}, newflags...; name="bands")

Uses the information from the template and supplied `kpoints` to generate a bands calculation.
Extra flags can be supplied which will be set for the generated calculation.
"""
function gencalc_bands(template::DFCalculation{QE}, kpoints::Vector{<:NTuple{4}},
                       newflags...; name = "bands")
    return calculation_from_kpoints(template, name, kpoints, :calculation => "bands",
                                    newflags...)
end

"""
    gencalc_nscf(template::DFCalculation{QE}, kpoints::NTuple{3, Int}, newflags...; name="nscf")

Uses the information from the template and supplied `kpoints` to generate an nscf calculation.
Extra flags can be supplied which will be set for the generated calculation.
"""
function gencalc_nscf(template::DFCalculation{QE}, kpoints::NTuple{3,Int}, newflags...;
                      name = "nscf")
    return calculation_from_kpoints(template, name, kpoints, :calculation => "nscf",
                                    newflags...)
end

"""
    gencalc_projwfc(template::DFCalculation{QE}, Emin, Emax, DeltaE, newflags...; name="projwfc")

Uses the information from the template and supplied `kpoints` to generate a `projwfc.x` calculation.
Extra flags can be supplied which will be set for the generated calculation.
"""
function gencalc_projwfc(template::DFCalculation{QE}, Emin, Emax, DeltaE, extraflags...;
                         name = "projwfc")
    occflag = flag(template, :occupations)
    ngauss  = 0
    if occflag == "smearing"
        smearingflag = flag(template, :smearing)
        if smearingflag ∈ ("methfessel-paxton", "m-p", "mp")
            ngauss = 1
        elseif smearingflag ∈ ("marzari-vanderbilt", "cold", "m-v", "mv")
            ngauss = -1
        elseif smearingflag ∈ ("fermi-dirac", "f-d", "fd")
            ngauss = -99
        end
    end
    tdegaussflag = flag(template, :degauss)
    degauss = tdegaussflag != nothing ? tdegaussflag : 0.0
    if length(execs(template)) == 2
        excs = [execs(template)[1],
                Exec(; exec = "projwfc.x", dir = execs(template)[end].dir)]
    else
        excs = [Exec(; exec = "projwfc.x", dir = execs(template)[end].dir)]
    end

    out = DFCalculation(deepcopy(template); name=name,  execs = excs, data=InputData[])
    set_name!(out, "projwfc")
    empty!(out.flags)
    set_flags!(out, :Emin => Emin, :Emax => Emax, :DeltaE => DeltaE, :ngauss => ngauss,
               :degauss => degauss; print = false)
    set_flags!(out, extraflags...)
    return out
end

"""
    gencalc_wan(nscf::DFCalculation{QE}, structure::AbstractStructure, Emin, wanflags...;
                Epad     = 5.0,
                wanexec  = Exec(exec="wannier90.x", dir=""))

Generates a Wannier90 calculation to follow on the supplied `nscf` calculation. It uses the projections defined in the `structure`, and starts counting the required amount of bands from `Emin`.
The `nscf` needs to have a valid output since it will be used in conjunction with `Emin` to find the required amount of bands and energy window for the Wannier90 calculation.
"""
function gencalc_wan(nscf::DFCalculation{QE}, structure::AbstractStructure, Emin,
                     wanflags...; Epad = 5.0,
                     wanexec = Exec(; exec = "wannier90.x", dir = ""))
    hasoutput_assert(nscf)
    iscalc_assert(nscf, "nscf")
    hasprojections_assert(structure)
    if iscolin(structure)
        wannames = ["wanup", "wandn"]
        @info "Spin polarized calculation found (inferred from nscf calculation)."
    else
        wannames = ["wan"]
    end

    if flag(nscf, :nosym) != true
        @info "'nosym' flag was not set in the nscf calculation.
                If this was not intended please set it and rerun the nscf calculation.
                This generally gives errors because of omitted kpoints, needed for pw2wannier90.x"
    end
    wanflags = wanflags != nothing ? SymAnyDict(wanflags) : SymAnyDict()
    if issoc(nscf)
        wanflags[:spinors] = true
    end

    nwann = nprojections(structure)
    @info "num_wann=$nwann (inferred from provided projections)."
    wanflags[:num_wann] = nwann
    kpoints = data(nscf, :k_points).data
    wanflags[:mp_grid] = kakbkc(kpoints)
    @info "mp_grid=$(join(wanflags[:mp_grid]," ")) (inferred from nscf calculation)."
    wanflags[:preprocess] = true
    isnoncolin(structure) && (wanflags[:spinors] = true)
    kdata = InputData(:kpoints, :none, [k[1:3] for k in kpoints])

    wancalculations = DFCalculation{Wannier90}[]
    for wanfil in wannames
        push!(wancalculations,
              DFCalculation{Wannier90}(name = wanfil, dir = dir(nscf), flags = copy(wanflags), data = [kdata],
                                       execs = [wanexec], run = true))
    end

    if length(wancalculations) > 1
        set_flags!(wancalculations[1], :spin => "up")
        set_flags!(wancalculations[2], :spin => "down")
    end

    map(x -> set_wanenergies!(x, structure, nscf, Emin; Epad = Epad), wancalculations)
    return wancalculations
end
kakbkc(kgrid) = length.(unique.([[n[i] for n in kgrid] for i in 1:3]))

"""
    gencalc_wan(nscf::DFCalculation{QE}, structure::AbstractStructure, projwfc::DFCalculation{QE}, threshold::Real, wanflags...; kwargs...)

Generates a wannier calculation, that follows on the `nscf` calculation. Instead of passing Emin manually, the output of a projwfc.x run
can be used together with a `threshold` to determine the minimum energy such that the contribution of the
projections to the DOS is above the `threshold`.
"""
function gencalc_wan(nscf::DFCalculation{QE}, structure::AbstractStructure,
                     projwfc::DFCalculation{QE}, threshold::Real, args...; kwargs...)
    hasexec_assert(projwfc, "projwfc.x")
    hasoutput_assert(projwfc)
    Emin = Emin_from_projwfc(structure, projwfc, threshold)
    return gencalc_wan(nscf, structure, Emin, args...; kwargs...)
end

"""
    isconverged(c::DFCalculation{QE})

Returns whether an `scf` calculation was converged.
"""
function isconverged(c::DFCalculation{QE})
    hasoutput_assert(c)
    iscalc_assert(c, "scf")
    return outputdata(c)[:converged]
end
