using DFControl, Base.Test
using Plots
# @test plot_qe_bands(joinpath(@__DIR__,"../assets/outputs/bands.out")) != nothing
# @test plot_qe_kpdos(joinpath(@__DIR__,"../assets/outputs/kpdos.out")) != nothing
test_bands = read_qe_bands_file(joinpath(@__DIR__,"../assets/outputs/bands.out"))
@test plot(test_bands[1],fermi=3)         != nothing
@test plot(test_bands[3],:relative_cart)  != nothing
@test plot(test_bands[4],:relative_cryst) != nothing
@test plot(test_bands)                    != nothing

t_fermi = read_fermi_from_qe_file(joinpath(@__DIR__,"../assets/outputs/scf.out"))
t_eigval = test_bands[1].eigvals[1]
@test apply_fermi_level(test_bands,3.2)[1].eigvals[1] == t_eigval - 3.2f0
@test apply_fermi_level(test_bands,joinpath(@__DIR__,"../assets/outputs/scf.out"))[1].eigvals[1] == t_eigval - t_fermi
apply_fermi_level!(test_bands,3.2)
@test test_bands[1].eigvals[1] == t_eigval - 3.2f0
t_eigval = test_bands[1].eigvals[1]
apply_fermi_level!(test_bands,joinpath(@__DIR__,"../assets/outputs/scf.out"))
@test test_bands[1].eigvals[1] == t_eigval - t_fermi
