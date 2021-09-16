using Colors
using .Plots
using LinearAlgebra
using RecipesBase
# using ..Client
# using ..Jobs
# using ..Structures

Base.:*(f::Number, r::RGB) = RGB(f * r.r, f * r.b, f * r.g)
Base.:+(r1::RGB, r2::RGB) = RGB(r1.r + r2.r, r1.b + r2.b, r1.g + r2.g)

function blend_color(contribs::Vector, at_colors)
    if any(isnan, contribs)
        return RGB(0.0, 0.0, 0.0)
    end
    if length(contribs) == 1
        return at_colors[1]
    end
    result = weighted_color_mean(normalize([contribs[1], contribs[2]])[1], at_colors[1],
                                 at_colors[2])
    totcontrib = contribs[1] + contribs[2]
    for i in 3:length(contribs)
        result = weighted_color_mean(normalize([totcontrib, contribs[i]])[1], result,
                                     at_colors[i])
        totcontrib += contribs[i]
    end
    return result
end

#@TODO create nice ticks when ks=:relative
# """
# Recipe to plot a Band. If ks is set to relative it calculates the relative offset of the k_points to the middle of the k-point path and puts that on the x-axis.

# Input:    band::Band
# Optional: ks=nothing
# Kwargs:   fermi=0, -> applies fermi level to band eigenvalues before plotting.
# linewidth=2
# """
@recipe function f(band::Band, ks = nothing; fermi = 0, linewidth = 2)
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
    title --> "Eigenvalues"
    yguide --> "Energy (eV)"
    legend --> false
    out = band.eigvals
    return ks, out
end

@recipe function f(bands::Vector{<:Band}, ks = nothing)
    for (i, band) in enumerate(bands)
        @series begin
            band, ks
        end
    end
end

@recipe function f(job::Job, ymin, ymax, occupy_ratio = 0.2; overlap_spin = false)
    palette_ = ismissing(Plots.default(:palette)) ? :default : Plots.default(:palette)
    tc = Plots.plot_color(pop!(plotattributes, :seriescolor,
                                              RGB.(Plots.color_list(Plots.palette(palette_)))))
    plt_colors = tc isa Colorant ? repeat([tc], 4) : repeat(tc, 4)
    ylims --> [ymin, ymax]
    gridalpha --> 0.9
    if !any(x -> eltype(x) == QE, job.calculations)
        error("output plotting only implemented for QE jobs.")
    end
    outdat = Client.outputdata(job)
    frmi = Jobs.readfermi(job, outdat)
    fermi --> frmi
    bands = Jobs.readbands(job, outdat)
    if bands === nothing
        error("No bands found in job $(job.name).")
    end

    # Bands part
    ks = Structures.high_symmetry_kpoints(job.structure)
    tick_vals = Int[]
    tick_syms = String[]
    kpoints = bands isa NamedTuple ? bands.up[1].k_points_cryst : bands[1].k_points_cryst
    for (i, k) in enumerate(kpoints)
        if ks !== nothing
            kpath = ks.kpoints
            for (sym, vec) in kpath
                if vec == k
                    push!(tick_vals, i)
                    push!(tick_syms, " " * string(sym) * " ")
                end
            end
        elseif norm(k) == 0
            push!(tick_vals, i)
            push!(tick_syms, " Γ ")
        end
    end
    if bands isa NamedTuple
        window_ids = map(x -> findall(x) do b
                             min = minimum(b.eigvals) - frmi
                             max = maximum(b.eigvals) - frmi
                             return !((min < ymin && max < ymin) ||
                                      (min > ymax && max > ymax))
                         end, bands)
    else
        window_ids = (findall(bands) do b
                          min = minimum(b.eigvals) - frmi
                          max = maximum(b.eigvals) - frmi
                          return !((min < ymin && max < ymin) || (min > ymax && max > ymax))
                      end,)
    end
    window_ids === nothing && error("No bands inside window")
    # We define a single band plotting series here
    function plot_band(band, color, label, subplot)
        if overlap_spin
            tit = ""
        else
            tit = subplot == 1 ? (length(bands) == 2 ? "Spin up" : "Eigenvalues") :
                  "Spin down"
        end
        @series begin
            xticks --> (tick_vals, tick_syms)
            title := tit
            yguide := subplot == 1 ? "Energy (eV)" : ""
            label := label
            subplot := subplot
            seriescolor := color
            legend := false
            1:length(kpoints), band.eigvals .- frmi
        end
    end

    # PDOS part
    projwfc = Utils.getfirst(x -> Calculations.isprojwfc(x) && haskey(outdat, x.name), job.calculations)
    if projwfc !== nothing
        if bands isa NamedTuple && !overlap_spin
            doswindow = 3
            layout --> (1, 3)
        else
            doswindow = 2
            layout --> (1, 2)
        end
        states, projbands = outdat[projwfc.name][:states], outdat[projwfc.name][:bands]
        # First we find the amount that all the states appear in the window
        state_occupations = zeros(length(states))
        for ib in 1:(bands isa NamedTuple ? 2 : 1)
            for wid in window_ids[ib]
                b = projbands[wid]
                for ψ in b.extra[:ψ]
                    for is in eachindex(states)
                        state_occupations[is] += ψ[is]
                    end
                end
            end
        end
        # Now we take the most occupied ones, and somehow find out where there's a sudden dropoff
        max_occ = maximum(state_occupations)
        # sorted_occ = sortperm(state_occupations, rev=true)
        # goodids = findall(i -> state_occupations[sorted_occ][i] > occupy_ratio * max_occ, 1:length(state_occupations))
        # ats_orbs = unique(map(x -> (job.structure.atoms[x.atom_id].name, orbital(x.l).name), states[sorted_occ][goodids]))
        goodids = findall(i -> state_occupations[i] > occupy_ratio * max_occ,
                          1:length(state_occupations))
        ats_orbs = unique(map(x -> (job.structure.atoms[x.atom_id].name, Structures.orbital(x.l)),
                              states[goodids]))
        @info "Found $(length(ats_orbs)) atomic orbitals that satisfy the minimum occupation:\n$ats_orbs"

        atom_colors = bands isa NamedTuple ?
                      [plt_colors[2:length(ats_orbs)+1],
                       plt_colors[length(ats_orbs)+2:2*length(ats_orbs)+1]] :
                      [plt_colors[2:length(ats_orbs)+1]]

        bands = bands isa NamedTuple ? bands : [bands]
        band_contribs = [[[zeros(length(ats_orbs)) for i in 1:length(kpoints)]
                          for i1 in 1:length(window_ids[d])] for d in 1:length(bands)]

        @info "Reading pdos files and generating band coloring..."
        for (ia, (atsym, orb)) in enumerate(ats_orbs)
            energies, pd =outdat[projwfc.name][:energies], outdat[projwfc.name][:pdos][atsym][orb]
            #Plots PDOS
            if size(pd, 2) == 2
                @series begin
                    label --> "$(atsym)_$(orb.name)_up"
                    yguide --> ""
                    subplot := doswindow
                    seriescolor := atom_colors[1][ia]
                    title := "DOS"
                    pd[:, 1], energies .- frmi
                end
                @series begin
                    label --> "$(atsym)_$(orb.name)_down"
                    yguide --> ""
                    subplot := doswindow
                    seriescolor := atom_colors[2][ia]
                    title := "DOS"
                    -1 .* pd[:, 2], energies .- frmi
                end
            else
                @series begin
                    label --> "$(atsym)_$(orb.name)"
                    yguide --> ""
                    subplot := doswindow
                    seriescolor := atom_colors[1][ia]
                    title := "DOS"
                    pd, energies .- frmi
                end
            end
            #Calculate band colors
            for (iud, (bnds, contribs)) in enumerate(zip(bands, band_contribs))
                for (ib, b) in enumerate(bnds[window_ids[iud]])
                    for ik in 1:length(kpoints)
                        ibin = findfirst(x -> energies[x] < b.eigvals[ik] <= energies[x+1],
                                         1:length(energies)-1)
                        contribs[ib][ik][ia] += ibin === nothing ? 0.0 : abs(pd[ibin, iud])
                    end
                end
            end
        end
        for contribs in band_contribs
            for ib in contribs
                for ik in 1:length(ib)
                    ib[ik] .= normalize(ib[ik])
                end
            end
        end
        band_colors = [[[blend_color(band_contribs[i][ib][ik], atom_colors[i])
                         for ik in 1:length(kpoints)] for ib in 1:length(window_ids[i])]
                       for i in 1:length(band_contribs)]
        @info "Plotting bands..."
        for (iplt, (bnds, colors)) in enumerate(zip(bands, band_colors))
            if length(bands) == 2
                lab = iplt == 1 ? "up" : "down"
            else
                lab = ""
            end
            for (ib, (b, c)) in enumerate(zip(bnds[window_ids[iplt]], colors))
                plot_band(b, c, ib == 1 ? lab : "", overlap_spin ? 1 : iplt)
            end
        end

        # If no pdos is present
    else
        if bands isa NamedTuple && !overlap_spin
            layout := (1, 2)
        else
            layout := (1, 1)
        end
        @info "Plotting bands..."
        if bands isa NamedTuple
            #loop over up down
            for (iplt, bnds) in enumerate(bands)
                lab = iplt == 1 ? "up" : "down"
                color = iplt == 1 ? :blue : :red
                #loop over bands inside window
                for (ib, b) in enumerate(bnds[window_ids[iplt]])
                    plot_band(b, color, ib == 1 ? lab : "", overlap_spin ? 1 : iplt)
                end
            end
        else
            for (ib, b) in enumerate(bands[window_ids[1]])
                plot_band(b, plt_colors[mod1(ib, length(plt_colors))], "", 1)
            end
        end
    end
end
