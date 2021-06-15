issoccalc(input::DFInput{Wannier90}) = flag(input, :spinors) == true

readoutput(input::DFInput{Wannier90}) = wan_read_output(outpath(input))

for f in (:cp, :mv)
    @eval function Base.$f(i::DFInput{Wannier90}, dest::String; kwargs...)
        for glob in ("$(name(i))", "UNK") # seedname should also cover generated pw2wannier90 files
            for f in searchdir(i, glob)
                $f(f, joinpath(dest, splitdir(f)[end]); kwargs...)
            end
        end
    end
end

"""
    set_kpoints!(input::DFInput, k_grid::NTuple{3, Int})

Sets the kpoints of the input.
"""
function set_kpoints!(input::DFInput{Wannier90}, k_grid::NTuple{3, Int}; print=true)
    set_flags!(input, :mp_grid => [k_grid...], print=print)
    set_data!(input, :kpoints, kgrid(k_grid..., :wan), print=print)
    return input
end

"""
    set_wanenergies!(waninput::DFInput{Wannier90}, structure::AbstractStructure, nscf::DFInput , Emin::Real; Epad=5.0)

Automatically calculates and sets the wannier energies. This uses the projections,
`Emin` and the output of the nscf calculation to infer the other limits.
`Epad` allows one to specify the padding around the inner and outer energy windows
"""
function set_wanenergies!(input::DFInput{Wannier90}, structure::AbstractStructure, nscf::DFInput, Emin::Real; Epad=5.0)
    hasoutput_assert(nscf)
    @assert isnscf(nscf) "Please provide a valid nscf calculation."
    hasprojections_assert(structure)

    bands = readbands(nscf)
    nwann = nprojections(structure)
    @info "num_wann=$nwann (inferred from provided projections)."

    if length(bands) == 2
        num_bands = length(bands[1])
        if hasflag(input, :spin) && input[:spin] == "up"
            winmin, frozmin, frozmax, winmax = wanenergyranges(Emin, nwann, bands.up, Epad)
        elseif hasflag(input, :spin) && input[:spin] == "down"
            winmin, frozmin, frozmax, winmax = wanenergyranges(Emin, nwann, bands.down, Epad)
        end
    else
        num_bands = length(bands)
        winmin, frozmin, frozmax, winmax = wanenergyranges(Emin, nwann, bands, Epad)
    end

    set_flags!(input, :dis_win_min => winmin, :dis_froz_min => frozmin, :dis_froz_max => frozmax, :dis_win_max => winmax, :num_wann => nwann, :num_bands=>num_bands;print=false)
    return input
end

