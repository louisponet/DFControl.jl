issoccalc(calculation::DFCalculation{Wannier90}) = flag(calculation, :spinors) == true

readoutput(calculation::DFCalculation{Wannier90}) = wan_read_output(outpath(calculation))

for f in (:cp, :mv)
    @eval function Base.$f(i::DFCalculation{Wannier90}, dest::String; kwargs...)
        for glob in ("$(name(i))", "UNK") # seedname should also cover generated pw2wannier90 files
            for f in searchdir(i, glob)
                $f(f, joinpath(dest, splitdir(f)[end]); kwargs...)
            end
        end
    end
end

"""
    set_kpoints!(calculation::DFCalculation, k_grid::NTuple{3, Int})

Sets the kpoints of the calculation.
"""
function set_kpoints!(calculation::DFCalculation{Wannier90}, k_grid::NTuple{3, Int}; print=true)
    set_flags!(calculation, :mp_grid => [k_grid...], print=print)
    set_data!(calculation, :kpoints, kgrid(k_grid..., :wan), print=print)
    return calculation
end

"""
    set_wanenergies!(wancalculation::DFCalculation{Wannier90}, structure::AbstractStructure, nscf::DFCalculation , Emin::Real; Epad=5.0)

Automatically calculates and sets the wannier energies. This uses the projections,
`Emin` and the output of the nscf calculation to infer the other limits.
`Epad` allows one to specify the padding around the inner and outer energy windows
"""
function set_wanenergies!(calculation::DFCalculation{Wannier90}, structure::AbstractStructure, nscf::DFCalculation, Emin::Real; Epad=5.0)
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
            winmin, frozmin, frozmax, winmax = wanenergyranges(Emin, nwann, bands.down, Epad)
        end
    else
        num_bands = length(bands)
        winmin, frozmin, frozmax, winmax = wanenergyranges(Emin, nwann, bands, Epad)
    end

    set_flags!(calculation, :dis_win_min => winmin, :dis_froz_min => frozmin, :dis_froz_max => frozmax, :dis_win_max => winmax, :num_wann => nwann, :num_bands=>num_bands;print=false)
    return calculation
end

