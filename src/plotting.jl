# not sure if we want to set defaults especially plot font
# plot_font=font(15,"DejaVu Sans")
# pyplot(lab="",yguidefont=plot_font,ytickfont=plot_font,xtickfont=plot_font,legendfont=plot_font)

#These should be in a plotting script.
# """
# Plots all the bands found in an output file of a Quantum Espresso bands calculation.
#
# Input:    filename::String,
# Optional: indices=nothing -> indices of the bands that you want to plot (low to high).
# """
# function plot_qe_bands(filename::String,indices=nothing,args...;kwargs...)
#   bands = read_qe_bands_file(filename)
#   if indices!=nothing
#     bands = bands[indices...]
#   end
#   plot(bands,args...;kwargs...)
# end
#
# #@Test: see if this is the correct way to use the tickvals!
# #@TODO: maybe put some good defaults in place?
# """
# Plots the k_resolved partial density of states from a Quantum Espresso projwfc output file.
# Only use this if the flag kresolveddos=true in the projwfc input file!!
# """
# function plot_qe_kpdos(filename::String, column=1 ; kwargs...)
#   array,(tickvals,ticks) = read_qe_kpdos(filename,column)
#   heatmap(array,yticks=(tickvals,ticks),kwargs...)
# end

#@TODO create nice ticks when ks=:relative
"""
Recipe to plot a DFBand. If ks is set to relative it calculates the relative offset of the k_points to the middle of the k-point path and puts that on the x-axis.

Input:    band::DFBand
Optional: ks=nothing
Kwargs:   fermi=0, -> applies fermi level to band eigenvalues before plotting.
          linewidth=2
"""
@recipe function f(band::DFBand, ks=nothing; fermi=0, linewidth=2)
  if ks==:relative_cart
    ks=[]
    k_m = band.k_points_cart[div(size(band.k_points_cart)[1]+1,2)]
    for k in band.k_points_cart
      push!(ks,norm(k-k_m))
    end
    ks[1:div(length(ks),2)]=-ks[1:div(length(ks),2)]
  elseif ks==:relative_cryst
    ks=[]
    k_m = band.k_points_cryst[div(size(band.k_points_cryst)[1]+1,2)]
    for k in band.k_points_cryst
      push!(ks,norm(k-k_m))
    end
    ks[1:div(length(ks),2)]=-ks[1:div(length(ks),2)]
  else
    ks = collect(1:length(band.k_points_cart))
  end
  if fermi != 0
    band = apply_fermi_level(band,fermi)
  end
  linewidth --> linewidth
  title --> "Eigenvalues"
  out = band.eigvals
  yguide -->(haskey(d,:yguide) ? d[:yguide] : "energy (eV)")
  ks,out
end

# @Test check if this still works!
"""
To plot multiple bands on one plot.
"""
@recipe function f(bands::Array{<:DFBand,1},ks=nothing)
  for (i,band) in enumerate(bands)
    @series begin
      label --> (haskey(d,:label) ? d[:label] : "Band $i")
      band,ks
    end
  end
end
