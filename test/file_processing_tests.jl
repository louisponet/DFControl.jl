using DFControl, Base.Test


scf_input     = read_qe_input(joinpath(@__DIR__,"../assets/inputs/qe/scf.in"),Float64)
bands_input   = read_qe_input(joinpath(@__DIR__,"../assets/inputs/qe/bands.in"),Float64)
projwfc_input = read_qe_input(joinpath(@__DIR__,"../assets/inputs/qe/projwfc.in"),Float64)

@test scf_input.control_blocks[:control][:calculation] == "'scf'"
@test scf_input.pseudos[:Te] == "Te.rel-pbesol-dn-kjpaw_psl.0.2.2.UPF"
@test scf_input.k_points[:option] == :automatic
@test scf_input.k_points[:nk3] == 10
@test scf_input.atoms[:Te] == Point3D{Float64}(0.523252856, 0.523252856, 0.523252856)

@test bands_input.control_blocks[:control][:calculation] == "'bands'"
@test bands_input.pseudos[:Te] == scf_input.pseudos[:Te]
@test bands_input.cell_param[:bohr] == scf_input.cell_param[:bohr]
@test bands_input.k_points[:option] == :crystal_b
@test bands_input.atoms[:Te] == scf_input.atoms[:Te]

@test projwfc_input.control_blocks[:projwfc][:kresolveddos]==true




 

