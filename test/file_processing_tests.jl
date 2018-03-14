using DFControl, Base.Test


scf_input, qestructure = read_qe_input(joinpath(@__DIR__,"../assets/inputs/qe/scf.in"));

bands_input   = read_qe_input(joinpath(@__DIR__,"../assets/inputs/qe/bands.in"))[1];
projwfc_input = read_qe_input(joinpath(@__DIR__,"../assets/inputs/qe/projwfc.in"), exec="projwfc.x"=>Dict{Symbol, Any}())[1];
display(projwfc_input);
@test get_flag(scf_input, :calculation) == "'scf'";
@test get_block(scf_input,:k_points).option == :automatic;
@test get_block(scf_input,:k_points).data[3] == 10;

@test get_block(bands_input, :k_points).option == :crystal_b;

@test get_flag(projwfc_input,:kresolveddos) == true

test_bands = read_qe_bands_file(joinpath(@__DIR__,"../assets/outputs/bands.out"));
@test length(test_bands) == 48
@test length(test_bands[1].k_points_cart) == 201
@test test_bands[1].k_points_cart == test_bands[2].k_points_cart;
@test test_bands[1].k_points_cart[3][2]==1.475098449380172
@test test_bands[1].k_points_cryst[5][1] == 0.5

test_ks_cart,test_ks_cryst = read_ks_from_qe_output(joinpath(@__DIR__,"../assets/outputs/bands.out"));
@test length(test_ks_cart) == 201
@test length(test_ks_cryst) == 201
@test test_ks_cryst == test_bands[1].k_points_cryst;
@test test_ks_cart  == test_bands[1].k_points_cart;

@test read_fermi_from_qe_output(joinpath(@__DIR__,"../assets/outputs/scf.out")) == 5.1371

wan_test, wanstructure = read_wannier_input(joinpath(@__DIR__,"../assets/inputs/wannier/wan.win"));
@test length(get_data(wan_test,:kpoints)) == prod(get_flag(wan_test,:mp_grid))
@test get_flag(wan_test,:write_rmn)

test_filename = joinpath(@__DIR__,"../assets/inputs/wannier/wan_test.win")
write_input(wan_test, wanstructure, test_filename);
@test wan_test.flags == read_wannier_input(test_filename)[1].flags
wan_test2 = read_wannier_input(test_filename)[1]
@test wan_test.data_blocks[1].data == wan_test2.data_blocks[1].data
@test print_flag(scf_input,:pseudo_dir) == print_flag(bands_input,:pseudo_dir)
remove_flags!(wan_test2, :dis_win_max, :dis_win_min)
@test get_flag(wan_test2,:dis_win_max) == get_flag(wan_test2,:dis_win_min)
rm(test_filename);

@test size(DFControl.gen_k_grid(10,10,10,:wan))[1] == size(DFControl.gen_k_grid(10,10,10,:nscf))[1]

if isdefined(:default_pseudo_dirs)
  pr_l = length(default_pseudo_dirs)
else
  pr_l =0
end
set_default_pseudo_dir(:default, "/test/test/test");
@test isdefined(:default_pseudo_dirs)
@test DFControl.get_default_pseudo_dirs()[:default] == "/test/test/test"
remove_default_pseudo_dir(:default);
@test length(keys(default_pseudo_dirs))==pr_l

pr_s = DFControl.get_default_server();

set_default_server("test/default");
@test DFControl.get_default_server() == "test/default"
if pr_s != ""
  set_default_server(pr_s)
end
if isdefined(:default_job_header)
  pr_h = default_job_header
end
set_default_job_header(["asdf","asdf"]);
set_default_job_header(["asdf","asdf"]);
@test default_job_header ==["asdf","asdf"]
if isdefined(:pr_h)
  set_default_job_header(pr_h)
end

set_default_input(scf_input, qestructure, :scf);
set_default_input(scf_input, qestructure, :scf);
set_default_input(bands_input, qestructure, :bands);

@test default_inputs[:scf][1].run_command == scf_input.run_command
remove_default_input(:scf)
remove_default_input(:bands)
