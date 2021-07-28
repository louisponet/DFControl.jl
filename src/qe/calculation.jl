function set_noncolin_flags!(i::DFCalculation{QE})
    return set_flags!(i, :npsin => 4, :noncolin => true; print = false)
end

function sanitize_flags!(c::DFCalculation{QE}, structure::AbstractStructure)
    if isvcrelax(c)
        #this is to make sure &ions and &cell are there in the calculation 
        !hasflag(c, :ion_dynamics) && set_flags!(c, :ion_dynamics => "bfgs"; print = false)
        !hasflag(c, :cell_dynamics) &&
            set_flags!(c, :cell_dynamics => "bfgs"; print = false)
    end
    #TODO add all the required flags
    if exec(c, "pw.x") !== nothing
        @assert hasflag(c, :calculation) "Please set the flag for calculation with name: $(name(c))"
    end
    # setting hubbard and magnetization flags
    set_hubbard_flags!(c, structure)
    set_starting_magnetization_flags!(c, structure)

    # setting hubbard flags 
    pseudo_dir = pseudo(atoms(structure)[1]).dir # Pseudos should be all sanitized by now
    set_flags!(c, :pseudo_dir => pseudo_dir; print = false)

    return convert_flags!(c)
end

isbands(c::DFCalculation{QE})   = flag(c, :calculation) == "bands"
isnscf(c::DFCalculation{QE})    = flag(c, :calculation) == "nscf"
isscf(c::DFCalculation{QE})     = flag(c, :calculation) == "scf"
isvcrelax(c::DFCalculation{QE}) = flag(c, :calculation) == "vc-relax"
isrelax(c::DFCalculation{QE})   = flag(c, :calculation) == "relax"

function ispw(c::DFCalculation{QE})
    return isbands(c) || isnscf(c) || isscf(c) || isvcrelax(c) || isrelax(c)
end

isprojwfc(c::DFCalculation{QE}) = hasexec(c, "projwfc.x")
ishp(c::DFCalculation{QE})      = hasexec(c, "hp.x")
issoc(c::DFCalculation{QE})     = flag(c, :lspinorb) == true

function ismagnetic(c::DFCalculation{QE})
    return (hasflag(c, :nspin) && c[:nspin] > 0.0) ||
           (hasflag(c, :total_magnetization) && c[:total_magnetization] != 0.0)
end

function readoutput(c::DFCalculation{QE}; kwargs...)
    return qe_read_output(c; kwargs...)
end

pseudodir(c::DFCalculation{QE}) = flag(c, :pseudo_dir)

function outfiles(c::DFCalculation{QE})
    files = [outpath(c)]
    for (is, fuzzies) in zip((isprojwfc, ishp), (("pdos",), ("Hubbard_parameters",)))
        if is(c)
            for f in fuzzies
                append!(files, searchdir(c, f))
            end
        end
    end
    return filter(ispath, unique(files))
end

function set_cutoffs!(c::DFCalculation{QE}, ecutwfc, ecutrho)
    return set_flags!(c, :ecutwfc => ecutwfc, :ecutrho => ecutrho)
end

function Emin_from_projwfc(structure::AbstractStructure, projwfc::DFCalculation{QE},
                           threshold::Number)
    hasoutput_assert(projwfc)
    hasexec_assert(projwfc, "projwfc.x")
    if !haskey(outdata(projwfc), :states)
        states, bands             = qe_read_projwfc(outpath(projwfc))
        outdata(projwfc)[:states] = states
        outdata(projwfc)[:bands]  = bands
    else
        states, bands = outdata(projwfc)[:states], outdata(projwfc)[:bands]
    end

    mask = zeros(length(states))
    for (atid, at) in enumerate(atoms(structure))
        projs = projections(at)

        if isempty(projs)
            continue
        end

        stateids = Int[]
        for proj in projs
            orb = orbital(proj)
            push!.((stateids,), findall(x -> x.atom_id == atid && x.l == orb.l, states))
        end
        mask[stateids] .= 1.0
    end
    Emin = -10000.0
    for b in bands
        ψ = mean(b.extra[:ψ])
        tot_relevant_occupation = dot(mask, ψ)

        if tot_relevant_occupation > threshold
            Emin = minimum(b.eigvals)
            break
        end
    end
    if Emin == -10000.0
        error("Couldn't find any band with occupation of relevant projections above $threshold, were any set in the structure?")
    end
    return Emin
end

#asserts
function iscalc_assert(i::DFCalculation{QE}, calc)
    @assert flag(i, :calculation) == calc "Please provide a valid '$calc' calculation."
end

function set_hubbard_flags!(c::DFCalculation{QE}, str::AbstractStructure{T}) where {T}
    u_ats = unique(atoms(str))
    isdftucalc = any(x -> dftu(x).U != 0 ||
                              dftu(x).J0 != 0.0 ||
                              sum(dftu(x).J) != 0 ||
                              sum(dftu(x).α) != 0, u_ats) || hasflag(c, :Hubbard_parameters)
    isnc = isnoncolin(str)
    if isdftucalc
        Jmap = map(x -> copy(dftu(x).J), u_ats)
        Jdim = maximum(length.(Jmap))
        Jarr = zeros(Jdim, length(u_ats))
        for (i, J) in enumerate(Jmap)
            diff = Jdim - length(J)
            if diff > 0
                for d in 1:diff
                    push!(J, zero(eltype(J)))
                end
            end
            Jarr[:, i] .= J
        end
        set_flags!(c, :lda_plus_u    => true, :Hubbard_U     => map(x -> dftu(x).U, u_ats),
                   :Hubbard_alpha => map(x -> dftu(x).α, u_ats),
                   :Hubbard_beta  => map(x -> dftu(x).β, u_ats), :Hubbard_J     => Jarr,
                   :Hubbard_J0    => map(x -> dftu(x).J0, u_ats); print = false)
        isnc && set_flags!(c, :lda_plus_u_kind => 1; print = false)
    else
        rm_flags!(c, :lda_plus_u, :lda_plus_u_kind, :Hubbard_U, :Hubbard_alpha,
                  :Hubbard_beta, :Hubbard_J, :Hubbard_J0, :U_projection_type; print = false)
    end
end

function set_starting_magnetization_flags!(c::DFCalculation{QE},
                                           str::AbstractStructure{T}) where {T}
    u_ats = unique(atoms(str))
    mags = magnetization.(u_ats)
    starts = T[]
    θs = T[]
    ϕs = T[]
    ismagcalc = ismagnetic(str)
    isnc = isnoncolin(str)
    if (ismagcalc && isnc) || (flag(c, :noncolin) !== nothing && flag(c, :noncolin))
        for m in mags
            tm = normalize(m)
            if norm(m) == 0
                push!.((starts, θs, ϕs), 0.0)
            else
                θ = acos(tm[3]) * 180 / π
                ϕ = atan(tm[2], tm[1]) * 180 / π
                start = norm(m)
                push!(θs, θ)
                push!(ϕs, ϕ)
                push!(starts, start)
            end
        end
        set_flags!(c, :noncolin => true; print = false)
        rm_flags!(c, :nspin; print = false)
    elseif ismagcalc
        for m in mags
            push!.((θs, ϕs), 0.0)
            if norm(m) == 0
                push!(starts, 0)
            else
                push!(starts, sign(sum(m)) * norm(m))
            end
        end
        set_flags!(c, :nspin => 2; print = false)
    end
    return set_flags!(c, :starting_magnetization => starts, :angle1 => θs, :angle2 => ϕs;
                      print = false)
end

for f in (:cp, :mv)
    @eval function Base.$f(i::DFCalculation{QE}, dest::String; kwargs...)
        $f(inpath(i), joinpath(dest, infilename(i)); kwargs...)
        if hasoutfile(i)
            $f(outpath(i), joinpath(dest, outfilename(i)); kwargs...)
        end
        if isprojwfc(i)
            for f in searchdir(i, "pdos")
                $f(f, joinpath(dest, splitdir(f)[end]); kwargs...)
            end
        elseif ishp(i)
            for f in searchdir(i, "Hubbard_parameters")
                $f(f, joinpath(dest, splitdir(f)[end]); kwargs...)
            end
        end
        #TODO add ph.x outfiles
    end
end

function pdos(c::DFCalculation{QE}, atsym::Symbol, magnetic::Bool, soc::Bool,
              filter_word = "")
    @assert isprojwfc(c) "Please specify a valid projwfc calculation."
    kresolved = hasflag(c, :kresolveddos) && calculation[:kresolveddos]
    files = filter(x -> occursin("($atsym)", x) &&
                            occursin("#", x) &&
                            occursin(filter_word, x), searchdir(dir(c), "pdos"))
    @assert !isempty(files) "No pdos files found in calculation directory $(dir(c))"
    files = joinpath.((c,), files)
    energies, = kresolved ? qe_read_kpdos(files[1]) : qe_read_pdos(files[1])
    atdos = magnetic && !soc ? zeros(size(energies, 1), 2) : zeros(size(energies, 1))
    if kresolved
        for f in files
            if magnetic && !occursin(".5", f)
                tu = qe_read_kpdos(f, 2)[2]
                td = qe_read_kpdos(f, 3)[2]
                atdos[:, 1] .+= reduce(+, tu; dims = 2) ./ size(tu, 2)
                atdos[:, 2] .+= reduce(+, td; dims = 2) ./ size(tu, 2)
                # elseif occursin(".5", f)
            else
                t = qe_read_kpdos(f, 1)[2]
                atdos .+= (reshape(reduce(+, t; dims = 2), size(atdos, 1)) ./ size(t, 2))
            end
        end
    else
        for f in files
            if magnetic && !occursin(".5", f)
                atdos .+= qe_read_pdos(f)[2][:, 1:2]
                # elseif occursin(".5", f)
            else
                atdos .+= qe_read_pdos(f)[2][:, 1]
            end
        end
    end
    return (energies = energies, pdos = atdos)
end

function set_kpoints!(c::DFCalculation{QE}, k_grid::NTuple{3,Int}; print = true) #nscf
    calc = flag(c, :calculation)
    print && calc != "nscf" && (@warn "Expected calculation to be 'nscf'.\nGot $calc.")
    set_data!(c, :k_points, kgrid(k_grid..., :nscf); option = :crystal, print = print)
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