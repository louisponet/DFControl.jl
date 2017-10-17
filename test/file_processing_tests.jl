using DFControl, Base.Test


scf_input     = read_qe_input(joinpath(@__DIR__,"../assets/inputs/qe/scf.in"),Float64)
bands_input   = read_qe_input(joinpath(@__DIR__,"../assets/inputs/qe/bands.in"),Float64)
projwfc_input = read_qe_input(joinpath(@__DIR__,"../assets/inputs/qe/projwfc.in"),Float64)

@test scf_input.control_blocks[:control][:calculation] == "'scf'"
@test scf_input.pseudos[:Te] == "Te.rel-pbesol-dn-kjpaw_psl.0.2.2.UPF"
@test scf_input.k_points[:option] == :automatic
@test scf_input.k_points[:nk3] == 10
@test scf_input.atoms[:Te] == [Point3D{Float64}(0.523252856, 0.523252856, 0.523252856)]

@test bands_input.control_blocks[:control][:calculation] == "'bands'"
@test bands_input.pseudos[:Te] == scf_input.pseudos[:Te]
@test bands_input.cell_param[:bohr] == scf_input.cell_param[:bohr]
@test bands_input.k_points[:option] == :crystal_b
@test bands_input.atoms[:Te] == scf_input.atoms[:Te]

@test projwfc_input.control_blocks[:projwfc][:kresolveddos]==true

test_bands = read_qe_bands_file(joinpath(@__DIR__,"../assets/outputs/bands.out"))
kpdos_test,(ticks,tickvals)  = read_qe_kpdos(joinpath(@__DIR__,"../assets/outputs/kpdos.out"))
@test length(test_bands) == 48
@test length(test_bands[1].k_points_cart) == 201
@test test_bands[1].k_points_cart == test_bands[2].k_points_cart
@test test_bands[1].k_points_cart[3][2]==1.4750984f0
@test test_bands[1].k_points_cryst[5][1] == 0.5f0
@test size(kpdos_test)   == (2002,201)
@test kpdos_test[302,3] == -0.17

test_ks_cart,test_ks_cryst = read_ks_from_qe_bands_file(joinpath(@__DIR__,"../assets/outputs/bands.out"))
@test length(test_ks_cart) == 201
@test length(test_ks_cryst) == 201
@test test_ks_cryst == test_bands[1].k_points_cryst
@test test_ks_cart  == test_bands[1].k_points_cart

@test read_fermi_from_qe_file(joinpath(@__DIR__,"../assets/outputs/scf.out")) == 5.1371f0

wan_test = read_wannier_input(joinpath(@__DIR__,"../assets/inputs/wannier/wan.win"))
@test length(wan_test.k_points) == prod(wan_test.control_blocks[:control][:mp_grid])
@test length(wan_test.atoms) == 3
@test wan_test.control_blocks[:control][:write_rmn]
@test wan_test.atoms[:Te] == [Point3D{Float32}(0.52325284,0.52325284,0.52325284)]

test_filename = joinpath(@__DIR__,"../assets/inputs/wannier/wan_test.win")
write_wannier_input(test_filename,wan_test)
@test wan_test.control_blocks[:control] == read_wannier_input(test_filename).control_blocks[:control]
wan_test2 = read_wannier_input(test_filename)
@test wan_test.k_points == wan_test2.k_points
@test wan_test.atoms == wan_test2.atoms
@test wan_test.cell_param == wan_test2.cell_param
rm(test_filename)
