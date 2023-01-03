include(joinpath(DEPS_DIR, "wannier90flags.jl"))
const WAN_FLAGS = _WAN_FLAGS()
issoc(c::Calculation{Wannier90}) = get(c, :spinors, false) 
flagtype(::Type{Wannier90}, flag) = haskey(WAN_FLAGS, flag) ? WAN_FLAGS[flag] : Nothing
flagtype(::Calculation{Wannier90}, flag) = flagtype(Wannier90, flag)

"""
    gencalc_wan(nscf::Calculation{QE}, structure::Structure, Emin, wanflags...;
                Epad     = 5.0,
                wanexec  = Exec(path="wannier90.x"))

Generates a Wannier90 calculation to follow on the supplied `nscf` calculation. It uses the projections defined in the `structure`, and starts counting the required amount of bands from `Emin`.
The `nscf` needs to have a valid output since it will be used in conjunction with `Emin` to find the required amount of bands and energy window for the Wannier90 calculation.
"""
function gencalc_wan(nscf::Calculation{QE}, structure::Structure, bands, Emin, wanflags...; Epad = 5.0,
                     wanexec = Exec(; name="wannier90", path = "wannier90.x"))
    projs = vcat(map(structure.atoms) do x
                     ps = x.projections
                     # @assert !isempty(ps) "Please first set projections for all atoms in the Structure."
                     return ps
                 end...)

    Structures.sanitize!(projs, Calculations.issoc(nscf))

    if Structures.iscolin(structure)
        wannames = ["wanup", "wandn"]
        @info "Spin polarized calculation found (inferred from nscf calculation)."
    else
        wannames = ["wan"]
    end

    wanflags = wanflags !== nothing ? Dict{Symbol,Any}(wanflags) : Dict{Symbol,Any}()
    if issoc(nscf)
        wanflags[:spinors] = true
    end

    nwann = issoc(nscf) ? 2*sum(length, projs) : sum(length, projs)
    @info "num_wann=$nwann (inferred from provided projections)."
    wanflags[:num_wann] = nwann
    kpoints = data(nscf, :k_points).data
    wanflags[:mp_grid] = length.(unique.([[n[i] for n in kpoints] for i in 1:3]))
    @info "mp_grid=$(join(wanflags[:mp_grid]," ")) (inferred from nscf calculation)."
    wanflags[:preprocess] = true
    Structures.isnoncolin(structure) && (wanflags[:spinors] = true)
    kdata = InputData(:kpoints, :none, [k[1:3] for k in kpoints])

    wancalculations = Calculation{Wannier90}[]
    for wanfil in wannames
        push!(wancalculations,
              Calculation(; name = wanfil, infile = wanfil * ".win", outfile = wanfil * ".wout",
                                     flags = copy(wanflags), data = [kdata],
                                     exec  = wanexec, run = true))
    end

    if length(wancalculations) > 1
        set_flags!(wancalculations[1], :spin => "up")
        set_flags!(wancalculations[2], :spin => "down")
    end

    map(x -> set_wanenergies!(x, structure, bands, Emin; Epad = Epad), wancalculations)
    return wancalculations
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
    d = data(calculation, :kpoints)
    if d !== nothing
        d.data = kgrid(k_grid..., calculation)
    else
        push!(calculation.data, InputData(:kpoints, :none, kgrid(k_grid..., calculation)))
    end
    return calculation
end

"""
    set_wanenergies!(wancalculation::Calculation{Wannier90}, structure::Structure, bands::Vector{Band}, Emin::Real; Epad=5.0)

Automatically calculates and sets the wannier energies. This uses the projections,
`Emin` and the output of the nscf calculation to infer the other limits.
`Epad` allows one to specify the padding around the inner and outer energy windows
"""
function set_wanenergies!(calculation::Calculation{Wannier90}, structure::Structure, bands,
                          Emin::Real; Epad = 5.0)
    nprojs = sum(x -> sum(length.(x.projections)), structure.atoms)
    nwann = Calculations.issoc(calculation) ? 2*nprojs : nprojs
    @assert nwann != 0 "Please specify projections first."
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

function Emin_from_projwfc(structure::Structure, states, bands::Vector{Band},
                           threshold::Number)
    mask = zeros(length(states))
    for (atid, at) in enumerate(structure.atoms)
        projs = at.projections

        if isempty(projs)
            continue
        end

        stateids = Int[]
        for proj in projs
            orb = proj.orbital
            push!.((stateids,), findall(x -> x.atom_id == atid && x.l == orb.l, states))
        end
        mask[stateids] .= 1.0
    end
    Emin = -10000.0
    for b in bands
        ψ = sum(b.extra[:ψ]) / length(b.extra[:ψ])
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

function remote_calcs(job, _calculation::Calculation{Wannier90})
    calcs = RemoteHPC.Calculation[]
    filename   = _calculation.infile
    should_run = _calculation.run
    nscf = getfirst(x -> Calculations.isnscf(x), job.calculations)
    
    @assert nscf !== nothing "No NSCF found to generate pw2wannier90 from."
    @assert eltype(nscf) == QE "Only QE based Wannier90 jobs are supported."

    pw2wan_exec = Exec(name = "", path=joinpath(dirname(nscf.exec), "pw2wannier90.x"), modules=nscf.exec.modules)

    preprocess   = get(_calculation, :preprocess, false)
    return [RemoteHPC.Calculation(_calculation.exec, "-pp $filename > $(_calculation.outfile)", preprocess || should_run),
            RemoteHPC.Calculation(pw2wan_exec, "-pd .true. < pw2wan_$(splitext(filename)[1]).in > pw2wan_$(splitext(filename)[1]).out", preprocess || should_run),
            RemoteHPC.Calculation(_calculation.exec, "$filename > $(_calculation.outfile)", should_run)]
end
