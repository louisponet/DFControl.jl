

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
    yguide    -->(haskey(d,:yguide) ? d[:yguide] : "energy (eV)")
    out = band.eigvals
    ks, out
end

"""
To plot multiple bands on one plot.
"""
@recipe function f(bands::Array{<:DFBand,1}, ks=nothing)
    for (i, band) in enumerate(bands)
        @series begin
            label --> (haskey(d, :label) ? d[:label] : "Band $i")
            band, ks
        end
    end
end
