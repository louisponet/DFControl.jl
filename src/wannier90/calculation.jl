function kgrid(na, nb, nc, ::Type{Wannier90})
    return reshape([(a, b, c)
                    for a in collect(range(0; stop = 1, length = na + 1))[1:end-1],
                        b in collect(range(0; stop = 1, length = nb + 1))[1:end-1],
                        c in collect(range(0; stop = 1, length = nc + 1))[1:end-1]],
                   (na * nb * nc))
end

function set_kpoints!(calculation::DFCalculation{Wannier90}, k_grid::NTuple{3,Int};
                      print = true)
    set_flags!(calculation, :mp_grid => [k_grid...]; print = print)
    set_data!(calculation, :kpoints, kgrid(k_grid..., calculation); print = print)
    return calculation
end

"""
    set_wanenergies!(wancalculation::DFCalculation{Wannier90}, structure::AbstractStructure, nscf::DFCalculation , Emin::Real; Epad=5.0)

Automatically calculates and sets the wannier energies. This uses the projections,
`Emin` and the output of the nscf calculation to infer the other limits.
`Epad` allows one to specify the padding around the inner and outer energy windows
"""
function set_wanenergies!(calculation::DFCalculation{Wannier90},
                          structure::AbstractStructure, nscf::DFCalculation, Emin::Real;
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

