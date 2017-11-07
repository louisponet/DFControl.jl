using DFControl, Base.Test


scf_input     = read_qe_input(joinpath(@__DIR__,"../assets/inputs/qe/scf.in"),Float64)
bands_input   = read_qe_input(joinpath(@__DIR__,"../assets/inputs/qe/bands.in"),Float64)
projwfc_input = read_qe_input(joinpath(@__DIR__,"../assets/inputs/qe/projwfc.in"),Float64)
display(projwfc_input)
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
@test display(test_bands[1]) == display([test_bands[1]])[1]
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

print_qe_flags(:electrons)
print_qe_namelists()
@test size(gen_k_grid(10,10,10,:wan))[1] == size(gen_k_grid(10,10,10,:nscf))[1]


add_default_pseudo_dir(:default,"/test/test/test")
@test isdefined(:default_pseudo_dirs)
@test DFWannier.get_default_pseudo_dirs()[:default] == "/test/test/test"
remove_default_pseudo_dir(:default)
@test length(keys(default_pseudo_dirs))==0

pr_s = DFWannier.get_default_server()

set_default_server("test/default")
@test DFWannier.get_default_server() == "test/default"
if pr_s != ""
  set_default_server(pr_s)
end
if isdefined(:default_job_header)
pr_h = default_job_header
end
set_default_job_header(["asdf","asdf"])
set_default_job_header(["asdf","asdf"])
@test default_job_header ==["asdf","asdf"] 
if pr_h != nothing
  set_default_job_header(pr_h)
end

set_default_input(scf_input,:scf)
set_default_input(scf_input,:scf)
set_default_input(bands_input,:bands)
@test default_inputs[:scf].run_command == scf_input.run_command 
remove_default_input(:scf)
remove_default_input(:bands)
