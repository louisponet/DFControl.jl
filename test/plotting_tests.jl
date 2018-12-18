using DFControl, Test
using Plots
test_bands = qe_read_bands_file(joinpath(@__DIR__, "testassets", "test_job", "nscf.out"))
@test plot(test_bands[1], fermi=3)         != nothing
@test plot(test_bands[3],:relative_cart)  != nothing
@test plot(test_bands[4],:relative_cryst) != nothing
@test plot(test_bands)                    != nothing

t_fermi = qe_read_fermi_from_output(joinpath(@__DIR__, "testassets", "test_job", "nscf.out"))
t_eigval = test_bands[1].eigvals[1]
@test DFControl.apply_fermi_level.(test_bands, 3.2)[1].eigvals[1] == t_eigval - 3.2
@test DFControl.apply_fermi_level.(test_bands, joinpath(@__DIR__, "testassets", "test_job", "nscf.out"))[1].eigvals[1] == t_eigval - t_fermi
DFControl.apply_fermi_level!.(test_bands, 3.2)
@test test_bands[1].eigvals[1] == t_eigval - 3.2
t_eigval = test_bands[1].eigvals[1]
DFControl.apply_fermi_level!.(test_bands, joinpath(@__DIR__, "testassets", "test_job", "nscf.out"))
@test test_bands[1].eigvals[1] == t_eigval - t_fermi
