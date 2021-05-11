using Colors

const PLOT_COLORS = [
    RGB(0.121569,0.466667,0.705882),
    RGB(1.000000,0.498039,0.054902),
    RGB(0.737255,0.741176,0.133333),
    RGB(0.580392,0.403922,0.741176),
    RGB(0.890196,0.466667,0.760784),
    RGB(0.498039,0.498039,0.498039),
    RGB(0.090196,0.745098,0.811765),
    RGB(0.839216,0.152941,0.156863),
    RGB(0.172549,0.627451,0.172549),
    RGB(0.227451,0.003922,0.513725),
    RGB(0.549020,0.337255,0.294118),
    RGB(0.000000,0.262745,0.003922),
    RGB(0.058824,1.000000,0.662745),
    RGB(0.368627,0.000000,0.250980),
    RGB(0.737255,0.737255,1.000000),
    RGB(0.847059,0.686275,0.635294),
    RGB(0.721569,0.000000,0.501961),
    RGB(0.000000,0.305882,0.325490),
    RGB(0.419608,0.396078,0.000000),
    RGB(0.490196,0.007843,0.000000)
]

Base.:*(f::Number, r::RGB) = RGB(f*r.r, f*r.b, f*r.g)
Base.:+(r1::RGB, r2::RGB) = RGB(r1.r + r2.r, r1.b + r2.b, r1.g+r2.g)


function blend_color(contribs::Vector)
    if any(isnan, contribs)
        return RGB(0.0,0.0,0.0)
    end
    if length(contribs) == 1
        return PLOT_COLORS[1]
    end
    result = weighted_color_mean(normalize([contribs[1], contribs[2]])[1], PLOT_COLORS[1], PLOT_COLORS[2])
    totcontrib = contribs[1] + contribs[2]
    for i = 3:length(contribs)
        result = weighted_color_mean(normalize([totcontrib, contribs[i]])[1], result, PLOT_COLORS[i])
        totcontrib += contribs[i]
    end
    return result
end

#@TODO create nice ticks when ks=:relative
"""
Recipe to plot a DFBand. If ks is set to relative it calculates the relative offset of the k_points to the middle of the k-point path and puts that on the x-axis.

Input:    band::DFBand
Optional: ks=nothing
Kwargs:   fermi=0, -> applies fermi level to band eigenvalues before plotting.
linewidth=2
"""
@recipe function f(band::DFBand, ks=nothing; fermi=0, linewidth=2)
    if ks == :relative_cart
        ks = []
        k_m = band.k_points_cart[div(size(band.k_points_cart)[1] + 1, 2)]
        for k in band.k_points_cart
            push!(ks, norm(k - k_m))
        end
        ks[1:div(length(ks), 2)] = -ks[1:div(length(ks), 2)]
    elseif ks == :relative_cryst
        ks = []
        k_m = band.k_points_cryst[div(size(band.k_points_cryst)[1] + 1, 2)]
        for k in band.k_points_cryst
            push!(ks, norm(k - k_m))
        end
        ks[1:div(length(ks), 2)] = -ks[1:div(length(ks), 2)]
    else
        ks = collect(1:length(band.k_points_cart))
    end
    if fermi != 0
        band = apply_fermi_level(band, fermi)
    end
    linewidth --> linewidth
    title     --> "Eigenvalues"
    yguide    --> "Energy (eV)"
    legend    --> false
    out = band.eigvals
    ks, out
end

"""
To plot multiple bands on one plot.
"""
@recipe function f(bands::Vector{<:DFBand}, ks=nothing)
    for (i, band) in enumerate(bands)
        @series begin
            band, ks
        end
    end
end

@recipe function f(job::DFJob, ymin, ymax, occupy_ratio=0.2; overlap_spin=false)
    ylims --> [ymin, ymax]
    if !isQEjob(job)
        error("output plotting only implemented for QE jobs.")
    end
    frmi = readfermi(job)
    fermi --> frmi
    bands = readbands(job)
    if bands === nothing
        error("No bands found in job $(job.name).")
    end

    # Bands part
    ks = high_symmetry_kpoints(job.structure)
    tick_vals = Int[]
    tick_syms = String[]
    kpoints = bands isa NamedTuple ? bands.up[1].k_points_cryst : bands[1].k_points_cryst
    for (i, k) in enumerate(kpoints)
        if ks!==nothing
            kpath = ks.kpoints
            for (sym, vec) in kpath
                if vec == k
                    push!(tick_vals, i)
                    push!(tick_syms, " " * string(sym) * " ")
                end
            end
        elseif norm(k) == 0
            push!(tick_vals, i)
            push!(tick_syms, " Gamma ")
        end
    end
    if bands isa NamedTuple
        window_ids = findall(bands.up) do b
            min = minimum(b.eigvals) - frmi
            max = maximum(b.eigvals) - frmi
            return !((min < ymin && max < ymin) || (min > ymax && max > ymax))
        end
    else
        window_ids = findall(bands) do b
            min = minimum(b.eigvals) - frmi
            max = maximum(b.eigvals) - frmi
            return !((min < ymin && max < ymin) || (min > ymax && max > ymax))
        end
    end
    window_ids === nothing && error("No bands inside window")

    # We define a single band plotting series here
    function plot_band(band, color, label, subplot)
        @series begin
            xticks --> (tick_vals, tick_syms)
            title --> "Eigenvalues"
            yguide := subplot == 1 ? "Energy (eV)" : "" 
            label := label 
            subplot := subplot
            seriescolor := color
            legend := true
            1:length(kpoints), band.eigvals .- frmi
        end
    end

    # PDOS part
    projwfc = getfirst(x -> isprojwfccalc(x) && hasoutfile(x), inputs(job))
    if projwfc !== nothing
        if bands isa NamedTuple && !overlap_spin
            doswindow = 3
            layout --> (1,3)
        else
            doswindow = 2
            layout --> (1,2)
        end 
        states, projbands = qe_read_projwfc(outpath(projwfc))
        # First we find the amount that all the states appear in the window
        state_occupations = zeros(length(states))
        for wid in window_ids
            b = projbands[wid]
            for ψ in b.extra[:ψ]
                for is in eachindex(states)
                    state_occupations[is] += ψ[is]
                end
            end
        end
        # Now we take the most occupied ones, and somehow find out where there's a sudden dropoff
        max_occ = maximum(state_occupations)
        sorted_occ = sortperm(state_occupations, rev=true)
        goodids = findall(i -> state_occupations[sorted_occ][i] > occupy_ratio * max_occ, 1:length(state_occupations))
        ats_orbs = unique(map(x -> (atoms(job)[x.atom_id].name, orbital(x.l).name), states[sorted_occ][goodids]))
        @info "Found $(length(ats_orbs)) atomic orbitals that satisfy the minimum occupation:\n$ats_orbs" 
        atom_colors = PLOT_COLORS[1:length(ats_orbs)]

        bands = bands isa NamedTuple ? bands : [bands]
        band_contribs = [[[zeros(length(ats_orbs)) for i = 1:length(kpoints)] for i1 = 1:length(window_ids)] for d = 1:length(bands)]

        @info "Reading pdos files and generating band coloring..."
        for (ia, (c, (atsym, orb))) in enumerate(zip(atom_colors, ats_orbs))
            energies, pd = pdos(job, atsym, "("*string(orb))
            
            #Plots PDOS
            if size(pd, 2) == 2
                @series begin
                    label --> "$(atsym)_$(orb)_up"
                    yguide --> ""
                    subplot := doswindow
                    seriescolor := 0.8 * c 
                    title := "DOS"
                    pd[:,1], energies .- frmi
                end
                @series begin
                    label --> "$(atsym)_$(orb)_down"
                    yguide --> ""
                    subplot := doswindow
                    seriescolor := c
                    title := "DOS"
                    pd[:,2], energies .- frmi
                end
            else
                @series begin
                    label --> "$(atsym)_$(orb)"
                    yguide --> ""
                    subplot := doswindow
                    seriescolor := c
                    title := "DOS"
                    pd, energies .- frmi
                end
            end

            #Calculate band colors
            for (iud, (bnds, contribs)) in enumerate(zip(bands, band_contribs))
                for (ib, b) in enumerate(bnds[window_ids])
                    for ik in 1:length(kpoints)
                        ibin = findfirst(x -> energies[x] < b.eigvals[ik] <= energies[x+1], 1:length(energies)-1)
                       
                        contribs[ib][ik][ia] += ibin === nothing ? 0.0 : pd[ibin, iud]
                    end
                end
            end
        end
        for contribs in band_contribs
            
             contribs .= [normalize.(contribs[ib]) for ib =1:length(window_ids)]
         end
        band_colors = [[[blend_color(band_contribs[i][ib][ik]) for ik = 1:length(kpoints)] for ib = 1:length(window_ids)] for i=1:length(band_contribs)]
        @info "Plotting bands..."
        for (iplt, (bnds, colors)) in enumerate(zip(bands, band_colors))
            if length(bands) == 2
                lab = iplt == 1 ? "up" : "down"
            else
                lab = ""
            end
            for (ib, (b, c)) in enumerate(zip(bnds[window_ids], colors))
                plot_band(b, c, ib == 1 ? lab : "", overlap_spin ? 1 : iplt)
            end
        end

    # If no pdos is present
    else
        if bands isa NamedTuple && !overlap_spin
            layout := (1,2)
        else
            layout := (1,1)
        end 
        @info "Plotting bands..."
        if bands isa NamedTuple
            #loop over up down
            for (iplt, bnds) in enumerate(bands)
                lab = iplt == 1 ? "up" : "down"
                color = iplt == 1 ? :blue : :red
                #loop over bands inside window
                for (ib, b) in enumerate(bnds[window_ids])
                    plot_band(b, color, ib == 1 ? lab : "", overlap_spin ? 1 : iplt)
                end
            end
        else
            for (ib, b) in enumerate(bands[window_ids])
                plot_band(b, PLOT_COLORS[mod1(ib,length(PLOT_COLORS))], "", 1)
            end
        end
    end
end
