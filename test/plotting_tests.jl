using DFControl, Test
using Plots
test_bands = qe_parse_bands_file(joinpath(@__DIR__, "testassets", "test_job", "nscf.out"))
@test plot(test_bands[1]; fermi = 3) != nothing
@test plot(test_bands[3], :relative_cart) != nothing
@test plot(test_bands[4], :relative_cryst) != nothing
@test plot(test_bands) != nothing

t_fermi = qe_parse_fermi_from_output(joinpath(@__DIR__, "testassets", "test_job",
                                             "nscf.out"))
