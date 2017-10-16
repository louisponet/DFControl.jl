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

 

