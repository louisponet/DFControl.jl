include(joinpath(DFC.DEPS_DIR, "wannier90flags.jl"))
const WAN_FLAGS = _WAN_FLAGS()
flagtype(::Type{Wannier90}, flag) = haskey(WAN_FLAGS, flag) ? WAN_FLAGS[flag] : Nothing
flagtype(::Calculation{Wannier90}, flag) = flagtype(Wannier90, flag)

"""
    gencalc_wan(nscf::Calculation{QE}, structure::Structure, Emin, wanflags...;
                Epad     = 5.0,
                wanexec  = Exec(exec="wannier90.x", dir=""))

Generates a Wannier90 calculation to follow on the supplied `nscf` calculation. It uses the projections defined in the `structure`, and starts counting the required amount of bands from `Emin`.
The `nscf` needs to have a valid output since it will be used in conjunction with `Emin` to find the required amount of bands and energy window for the Wannier90 calculation.
"""
function gencalc_wan(nscf::Calculation{QE}, structure::Structure, Emin,
                     wanflags...; Epad = 5.0,
                     wanexec = Exec(; exec = "wannier90.x", dir = ""))
    hasoutput_assert(nscf)
    @assert nscf[:calculation] == "nscf"
    projs = vcat(map(structure.atoms) do x
        ps = x.projections
        @assert !isempty(ps) "Please first set projections for all atoms in the Structure."
        return ps
    end
    )

    DFC.Structures.sanitize!(projs, DFC.Calculation.issoc(nscf))
    
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
    wanflags = wanflags != nothing ? Dict{Symbol,Any}(wanflags) : Dict{Symbol,Any}()
    if issoc(nscf)
        wanflags[:spinors] = true
    end

    nwann = sum(length, projs)
    @info "num_wann=$nwann (inferred from provided projections)."
    wanflags[:num_wann] = nwann
    kpoints = data(nscf, :k_points).data
    wanflags[:mp_grid] = length.(unique.([[n[i] for n in kgrid] for i in 1:3]))
    @info "mp_grid=$(join(wanflags[:mp_grid]," ")) (inferred from nscf calculation)."
    wanflags[:preprocess] = true
    isnoncolin(structure) && (wanflags[:spinors] = true)
    kdata = InputData(:kpoints, :none, [k[1:3] for k in kpoints])

    wancalculations = Calculation{Wannier90}[]
    for wanfil in wannames
        push!(wancalculations,
              Calculation{Wannier90}(name = wanfil, dir = nscf.dir, flags = copy(wanflags), data = [kdata],
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
    gencalc_wan(nscf::Calculation{QE}, structure::Structure, projwfc::Calculation{QE}, threshold::Real, wanflags...; kwargs...)

Generates a wannier calculation, that follows on the `nscf` calculation. Instead of passing Emin manually, the output of a projwfc.x run
can be used together with a `threshold` to determine the minimum energy such that the contribution of the
projections to the DOS is above the `threshold`.
"""
function gencalc_wan(nscf::Calculation{QE}, structure::Structure,
                     projwfc::Calculation{QE}, threshold::Real, args...; kwargs...)
    hasexec_assert(projwfc, "projwfc.x")
    hasoutput_assert(projwfc)
    Emin = Emin_from_projwfc(structure, projwfc, threshold)
    return gencalc_wan(nscf, structure, Emin, args...; kwargs...)
end

"""
    gencalc_wan(job::Job, min_window_determinator::Real, extra_wan_flags...; kwargs...)

Automates the generation of wannier calculations based on the `job`.
When a projwfc calculation is present in the `job`, `min_window_determinator` will be used to
determine the threshold value for including a band in the window based on the projections, otherwise
it will be used as the `Emin` value from which to start counting the number of bands needed for all
projections.
`extra_wan_flags` can be any extra flags for the Wannier90 calculation such as `write_hr` etc.
"""
function gencalc_wan(job::Job, min_window_determinator::Real, extra_wan_flags...;
                     kwargs...)
    nscf_calculation = getfirst(x -> isnscf(x), job.calculations)
    projwfc_calculation = getfirst(x -> isprojwfc(x), job.calculations)
    if projwfc_calculation === nothing || !hasoutput(projwfc_calculation)
        @info "No projwfc calculation found with valid output, using $min_window_determinator as Emin"
        return gencalc_wan(nscf_calculation, job.structure, min_window_determinator,
                           extra_wan_flags...; kwargs...)
    else
        @info "Valid projwfc output found, using $min_window_determinator as the dos threshold."
        return gencalc_wan(nscf_calculation, job.structure, projwfc_calculation,
                           min_window_determinator, extra_wan_flags...; kwargs...)
    end
end

## Wannier90
function kgrid(na, nb, nc, ::Type{Wannier90})
    return reshape([(a, b, c)
                    for a in collect(range(0; stop = 1, length = na + 1))[1:end-1],
                        b in collect(range(0; stop = 1, length = nb + 1))[1:end-1],
                        c in collect(range(0; stop = 1, length = nc + 1))[1:end-1]],
                   (na * nb * nc))
end

function set_kpoints!(calculation::Calculation{Wannier90}, k_grid::NTuple{3,Int};
                      print = true)
    set_flags!(calculation, :mp_grid => [k_grid...]; print = print)
    set_data!(calculation, :kpoints, kgrid(k_grid..., calculation); print = print)
    return calculation
end

"""
    set_wanenergies!(wancalculation::Calculation{Wannier90}, structure::Structure, nscf::Calculation , Emin::Real; Epad=5.0)

Automatically calculates and sets the wannier energies. This uses the projections,
`Emin` and the output of the nscf calculation to infer the other limits.
`Epad` allows one to specify the padding around the inner and outer energy windows
"""
function set_wanenergies!(calculation::Calculation{Wannier90},
                          structure::Structure, nscf::Calculation, Emin::Real;
                          Epad = 5.0)
    hasoutput_assert(nscf)
    @assert isnscf(nscf) "Please provide a valid nscf calculation."
    hasprojections_assert(structure)

    bands = readbands(nscf)
    nwann = nprojections(structure)
    @info "num_wann=$nwann (inferred from provided projections)."

    if length(bands) == 2
        num_bands = length(bands[1])
        if hasflag(calculation, :spin) && calculation[:spin] == "up"
            winmin, frozmin, frozmax, winmax = wanenergyranges(Emin, nwann, bands.up, Epad)
        elseif hasflag(calculation, :spin) && calculation[:spin] == "down"
            winmin, frozmin, frozmax, winmax = wanenergyranges(Emin, nwann, bands.down,
                                                               Epad)
        end
    else
        num_bands = length(bands)
        winmin, frozmin, frozmax, winmax = wanenergyranges(Emin, nwann, bands, Epad)
    end

    set_flags!(calculation, :dis_win_min => winmin, :dis_froz_min => frozmin,
               :dis_froz_max => frozmax, :dis_win_max => winmax, :num_wann => nwann,
               :num_bands => num_bands; print = false)
    return calculation
end

#TODO: there is probably a slightly more optimal way for the frozen window, by checking max and min
#      such that every k-point has nbands inbetween.
function Emax(Emin, nbnd, bands)
    nbndfound = 0
    max = 0
    for b in bands
        if maximum(b.eigvals) >= Emin && nbndfound <= nbnd
            nbndfound += 1
            #maximum of allowed frozen window is the minimum of the first band>nbnd
            max = minimum(b.eigvals) - 0.005
        end
    end

    nbndfound < nbnd &&
        error("Number of needed bands for the projections ($nbnd) exceeds the amount of bands starting from \nEmin=$Emin ($nbndfound).\nRerun nscf with nbnd=$(length(bands) + nbnd - nbndfound).")
    return max
end

function wanenergyranges(Emin, nbnd, bands, Epad = 5)
    max = Emax(Emin, nbnd, bands)
    return (Emin - Epad, Emin, max, max + Epad)
end

function Emin_from_projwfc(structure::DFC.Structure, projwfc::Calculation{QE},
                           threshold::Number)
    DFC.hasoutput_assert(projwfc)
    @assert any(x -> x.exec == "projwfc.x", projwfc.execs)
    if !haskey(projwfc.outdata, :states)
        states, bands             = qe_read_projwfc(outpath(projwfc))
        projwfc.outdata[:states] = states
        projwfc.outdata[:bands]  = bands
    else
        states, bands = projwfc.outdata[:states], projwfc.outdata[:bands]
    end

    mask = zeros(length(states))
    for (atid, at) in enumerate(structure.atoms)
        projs = at.projections

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

