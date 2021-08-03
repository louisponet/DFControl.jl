#----------- Basic Interaction --------------#
"""
    data(calculation::Calculation, n::Symbol)

The former returns `calculation.data`, the later -- the `InputData` with name `n`.
"""
data(calculation::Calculation, n::Symbol) = getfirst(x -> x.name == n, calculation.data)

"""
    set_data!(calculation::Calculation, data::InputData)

If an `InputData` with the same name as `data` is already in `calculation`, it will be overwritten. Otherwise `data` gets pushed to the list of `InputData` blocks.
"""
function set_data!(calculation::Calculation, data::InputData)
    id = findfirst(x -> x.name == data.name, calculation.data)
    if id === nothing
        push!(calculation.data, data)
    else
        calculation.data[id] = data
    end
    return calculation
end

"""
    set_data!(calculation::Calculation, block_name::Symbol, new_block_data; option::Symbol=nothing, print::Bool=true)

Searches for an `InputData` for which `InputData.name == block_name`, and sets `Calculation.data = new_block_data`.
If `option` is specified it is set, i.e. `InputData.option = option`.
"""
function set_data!(calculation::Calculation, block_name::Symbol, new_block_data;
                   option::Symbol = nothing, print::Bool = true)
    setd = false
    for data_block in calculation.data
        if data_block.name == block_name
            if typeof(data_block.data) != typeof(new_block_data)
                if print
                    @warn "Overwritten data of type '$(typeof(data_block.data))' with type '$(typeof(new_block_data))'."
                end
            end
            old_data = data_block.data
            data_block.data = new_block_data
            data_block.option = option == nothing ? data_block.option : option
            if print
                @info "Block data '$(data_block.name)' in calculation  '$(calculation.name)' is now:\n\t$(string(data_block.data)) \n\toption: $(data_block.option)\n"
            end
            setd = true
        end
    end
    if !setd
        set_data!(calculation, InputData(block_name, option, new_block_data))
        setd = true
    end
    return setd, calculation
end

#---------- Extended Interaction ------------#

"""
    readbands(calculation::Calculation)

Parses the `outputdata` associated with `calculation` and returns the bands if present.

    readbands(job::Job)

Will try to find the first bandstructure calculation present in the `job` with a valid
output file, and return the `bands` if present.
If no normal bandstructure calculation is present, it will return the bands produced
by a possible `nscf` calculation if one exists with a valid output file.
"""
function readbands(calculation::Calculation)
    hasoutput_assert(calculation)
    to = readoutput(calculation)
    if haskey(to, :bands)
        return to[:bands]
    elseif haskey(to, :bands_up)
        return (up = to[:bands_up], down = to[:bands_down])
    else
        error("No bands found in $(calculation.name).")
    end
end

"""
    readfermi(calculation::Calculation)

Parses the `outputdata` associated with `calculation` and returns the fermi level if present.

    readfermi(job::Job)

Finds the first `scf` calculation present in the `job`,
reads its `outputdata` and returns the fermi level if present.
"""
function readfermi(calculation::Calculation)
    hasoutput_assert(calculation)
    to = readoutput(calculation)
    if haskey(to, :fermi)
        return to[:fermi]
    else
        error("No fermi found in $(calculation.name).")
    end
end

"""
    set_kpoints!(calculation::Calculation{QE}, k_grid::NTuple{3, Int}; print=true)
    set_kpoints!(calculation::Calculation{QE}, k_grid::NTuple{6, Int}; print=true)
    set_kpoints!(calculation::Calculation{QE}, k_grid::Vector{<:NTuple{4}}; print=true, k_option=:crystal_b)

Convenience function to set the `:k_points` data block of `calculation`.
The three different methods are targeted at `nscf`, `scf` or `vcrelax`,
and `bands` calculations, respectively.
For the `nscf` version an explicit list of `k_points` will be generated.

    set_kpoints!(calculation::Calculation{Wannier90}, k_grid::NTuple{3, Int})

Similar to the `nscf` targeted function in the sense that it will generate
an explicit list of `k_points`, adhering to the same rules as for the `nscf`.
The `mp_grid` flag will also automatically be set.
"""
function set_kpoints!(::Calculation{P}, args...; kwargs...) where {P}
    @error "set_kpoints! not implemented for package $P."
end

#-------- Generating new Calculations ---------- #

function calculation_from_kpoints(template::Calculation, newname, kpoints, newflags...)
    newcalc = Calculation(deepcopy(template); name = newname)
    set_flags!(newcalc, newflags...; print=false)
    set_name!(newcalc, newname)
    set_kpoints!(newcalc, kpoints; print = false)
    return newcalc
end

function gencalc_nscf(::Calculation{P}, args...) where {P}
    @error "gencalc_nscf is not implemented for package $P."
end

function gencalc_scf(::Calculation{P}, args...) where {P}
    @error "gencalc_scf is not implemented for package $P."
end

function gencalc_projwfc(::Calculation{P}, args...) where {P}
    @error "gencalc_projwfc is not implemented for package $P."
end

function gencalc_wan(::Calculation{P}, args...) where {P}
    @error "gencalc_wan is not implemented for package $P."
end

function gencalc_vcrelax(::Calculation{P}, args...) where {P}
    @error "gencalc_vcrelax is not implemented for package $P."
end

function gencalc_bands(::Calculation{P}, args...) where {P}
    @error "gencalc_bands is not implemented for package $P."
end

function isconverged(::Calculation{P}) where {P}
    @error "isconverged is not implemented for package $P."
end

"""
    kgrid(na, nb, nc, calculation)

Returns an array of k-grid points that are equally spaced, calculation can be either `:wan` or `:nscf`, the returned grids are appropriate as calculations for wannier90 or an nscf calculation respectively.
"""
kgrid(na, nb, nc, ::Calculation{T}) where {T} = kgrid(na, nb, nc, T)

## QE
function kgrid(na, nb, nc, ::Type{QE})
    return reshape([(a, b, c, 1 / (na * nb * nc))
                    for a in collect(range(0; stop = 1, length = na + 1))[1:end-1],
                        b in collect(range(0; stop = 1, length = nb + 1))[1:end-1],
                        c in collect(range(0; stop = 1, length = nc + 1))[1:end-1]],
                   (na * nb * nc))
end

function set_kpoints!(c::Calculation{QE}, k_grid::NTuple{3,Int}; print = true) #nscf
    print && !DFC.isnscf(c) && (@warn "Expected calculation to be 'nscf'.\nGot $calc.")
    set_data!(c, :k_points, kgrid(k_grid..., c); option = :crystal, print = print)
    prod(k_grid) > 100 && set_flags!(c, :verbosity => "high"; print = print)
    return c
end

function set_kpoints!(c::Calculation{QE}, k_grid::NTuple{6,Int}; print = true) #scf
    print && !(DFC.isscf(c) || DFC.isvcrelax(c) || DFC.isrelax(c)) && (@warn "Expected calculation to be scf, vc-relax, relax.\nGot $calc.")
    set_data!(c, :k_points, [k_grid...]; option = :automatic, print = print)
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
    set_data!(c, :k_points, k_grid; option = k_option, print = print)
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
function gencalc_bands(template::Calculation{QE}, kpoints::Vector{<:NTuple{4}},
                       newflags...; name = "bands")
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

    out = Calculation(deepcopy(template); name=name,  execs = excs, data=InputData[])
    set_name!(out, "projwfc")
    empty!(out.flags)
    set_flags!(out, :Emin => Emin, :Emax => Emax, :DeltaE => DeltaE, :ngauss => ngauss,
               :degauss => degauss; print = false)
    set_flags!(out, extraflags...)
    return out
end

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
    wanflags = wanflags != nothing ? SymAnyDict(wanflags) : SymAnyDict()
    if issoc(nscf)
        wanflags[:spinors] = true
    end

    nwann = sum(length, projs)
    @info "num_wann=$nwann (inferred from provided projections)."
    wanflags[:num_wann] = nwann
    kpoints = data(nscf, :k_points).data
    wanflags[:mp_grid] = kakbkc(kpoints)
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
kakbkc(kgrid) = length.(unique.([[n[i] for n in kgrid] for i in 1:3]))

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
    DFC.hasexec_assert(projwfc, "projwfc.x")
    if !haskey(projwfc.outdata, :states)
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

