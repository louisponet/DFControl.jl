using DFControl, Base.Test


scf_input   = read_qe_input(joinpath(@__DIR__,"/../assets/inputs/qe/scf.in")
bands_input = read_qe_input(joinpath(@__DIR__,"/../assets/inputs/qe/bands.in")
@test scf_input.control_blocks[:control][:calculation] == "'scf'"
@test bands_input.control_blocks[:control][:calculation] == "'scf'"
 

