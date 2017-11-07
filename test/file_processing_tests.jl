using DFControl, Base.Test


scf_input     = read_qe_input(joinpath(@__DIR__,"../assets/inputs/qe/scf.in"),Float64)
bands_input   = read_qe_input(joinpath(@__DIR__,"../assets/inputs/qe/bands.in"),Float64)
projwfc_input = read_qe_input(joinpath(@__DIR__,"../assets/inputs/qe/projwfc.in"),Float64)

@test get_flag(scf_input,:calculation) == "'scf'"
@test get_data(scf_input,:atomic_species)[:Te] == "Te.rel-pbesol-dn-kjpaw_psl.0.2.2.UPF"
@test get_block(scf_input,:k_points).option == :automatic
@test get_block(scf_input,:k_points).data[3] == 10
@test get_block(scf_input,:atomic_positions).data[:Te] == [Point3D{Float64}(0.523252856, 0.523252856, 0.523252856)]

@test get_data(scf_input,:atomic_species)[:Te] == get_data(bands_input,:atomic_species)[:Te]
@test get_block(bands_input, :cell_parameters).data == get_block(scf_input, :cell_parameters).data
@test get_block(bands_input, :k_points).option == :crystal_b
@test get_block(bands_input,:atoms_positions) == get_block(scf_input,:atoms_positions)

@test get_flag(projwfc_input,:kresolveddos) == true

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
@test length(get_data(wan_test,:kpoints)) == prod(get_flag(wan_test,:mp_grid))
@test length(get_data(wan_test,:atoms_frac)) == 2
@test get_flag(wan_test,:write_rmn)
@test get_data(wan_test,:atoms_frac)[:Te] == [Point3D{Float32}(0.52325284,0.52325284,0.52325284)]

test_filename = joinpath(@__DIR__,"../assets/inputs/wannier/wan_test.win")
write_input(wan_test,test_filename)
@test wan_test.flags == read_wannier_input(test_filename).flags
wan_test2 = read_wannier_input(test_filename)
@test wan_test.data_blocks[1].data == wan_test2.data_blocks[1].data
@test print_flag(scf_input,:pseudo_dir) == print_flag(bands_input,:pseudo_dir)
remove_flags!(wan_test2,[:dis_win_max,:dis_win_min])
@test get_flag(wan_test2,:dis_win_max)==get_flag(wan_test2,:dis_win_min)
rm(test_filename)
