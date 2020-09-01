

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

@recipe function f(job::DFJob, ymin, ymax, occupy_ratio=0.2)
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

    if bands isa NamedTuple
        layout := (1,2)
    else
        layout := (1,2)
    end 
    # Bands part
    kpath = high_symmetry_kpoints(job.structure)
    tick_vals = Int[]
    tick_syms = String[]
    kpoints = bands isa NamedTuple ? bands.up[1].k_points_cryst : bands[1].k_points_cryst
    for (i, k) in enumerate(kpoints)
        for (sym, vec) in kpath
            if vec == k
                push!(tick_vals, i)
                push!(tick_syms, string(sym))
            end
        end
    end
    subplot_counter = 1
    if bands isa NamedTuple
        @series begin
            ticks --> (tick_vals, tick_syms)
            label := "up"
            subplot := subplot_counter
            seriescolor := :red
            legend := true
            bands.up
        end
        @series begin
            ticks --> (tick_vals, tick_syms)
            label := "down"
            subplot := subplot_counter
            seriescolor := :blue
            marker := :dot
            legend := true
            bands.down
        end
        subplot_counter += 1
    else
        @series begin
            ticks --> (tick_vals, tick_syms)
            title := "Eigenvalues"
            subplot := 1
            bands
        end
        subplot_counter += 1
    end
    @show subplot_counter
    # PDOS part
    projwfc = getfirst(x -> isprojwfccalc(x) && hasoutfile(x), inputs(job))
    states, projbands = qe_read_projwfc(outpath(projwfc))
    window_ids = bands isa NamedTuple ? findall(x -> minimum(x.eigvals .- frmi) > ymin && maximum(x.eigvals .- frmi) < ymax, bands.up) : findall(x -> minimum(x.eigvals .- frmi) > ymin && maximum(x.eigvals .- frmi) < ymax, bands)
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
    
    for (atsym, orb) in ats_orbs
        energies, pd = pdos(job, atsym, string(orb))
        if size(pd, 2) == 2
            @series begin
                label --> "$(atsym)_$(orb)_up"
                subplot := subplot_counter
                pd[:,1], energies .- frmi
            end
            @series begin
                label --> "$(atsym)_$(orb)_down"
                subplot := subplot_counter
                pd[:,2], energies .- frmi
            end
        else
            @series begin
                label --> "$(atsym)_$(orb)"
                subplot := subplot_counter
                pd, energies .- frmi
            end
        end
    end
end
