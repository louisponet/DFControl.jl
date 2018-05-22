using DFControl, Base.Test


scf_input, qestructure = read_qe_input(joinpath(@__DIR__,"../assets/inputs/qe/scf.in"));

bands_input   = read_qe_input(joinpath(@__DIR__,"../assets/inputs/qe/bands.in"))[1];
projwfc_input = read_qe_input(joinpath(@__DIR__,"../assets/inputs/qe/projwfc.in"), exec=Exec("projwfc.x","",Dict{Symbol, Any}()))[1];
display(projwfc_input);
@test flag(scf_input, :calculation) == "'scf'";
@test block(scf_input,:k_points).option == :automatic;
@test block(scf_input,:k_points).data[3] == 10;

@test block(bands_input, :k_points).option == :crystal_b;

@test flag(projwfc_input,:kresolveddos) == true

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
@test length(data(wan_test,:kpoints)) == prod(flag(wan_test,:mp_grid))
@test flag(wan_test,:write_rmn)

test_filename = joinpath(@__DIR__,"../assets/inputs/wannier/wan_test.win")
save(wan_test, wanstructure, test_filename);
@test wan_test.flags == read_wannier_input(test_filename)[1].flags
wan_test2 = read_wannier_input(test_filename)[1]
@test wan_test.data[1].data == wan_test2.data[1].data
@test print_flag(scf_input,:pseudo_dir) == print_flag(bands_input,:pseudo_dir)
rmflags!(wan_test2, :dis_win_max, :dis_win_min)
@test flag(wan_test2,:dis_win_max) == flag(wan_test2,:dis_win_min)
rm(test_filename);

@test size(DFControl.kgrid(10,10,10,:wan))[1] == size(DFControl.kgrid(10,10,10,:nscf))[1]

if isdefined(:default_pseudo_dirs)
  pr_l = length(default_pseudo_dirs)
else
  pr_l =0
end
setdefault_pseudodir(:default, "/test/test/test");
@test isdefined(:default_pseudo_dirs)
@test DFControl.getdefault_pseudodirs()[:default] == "/test/test/test"
removedefault_pseudodir(:default);
@test length(keys(default_pseudo_dirs))==pr_l

pr_s = DFControl.getdefault_server();

setdefault_server("test/default");
@test DFControl.getdefault_server() == "test/default"
if pr_s != ""
  setdefault_server(pr_s)
end
if isdefined(:default_job_header)
  pr_h = default_job_header
end
setdefault_jobheader(["asdf","asdf"]);
setdefault_jobheader(["asdf","asdf"]);
@test default_job_header ==["asdf","asdf"]
if isdefined(:pr_h)
  setdefault_jobheader(pr_h)
end

setdefault_input(scf_input, qestructure, :scf);
setdefault_input(scf_input, qestructure, :scf);
setdefault_input(bands_input, qestructure, :bands);

@test default_inputs[:scf][1].runcommand.exec == scf_input.runcommand.exec
removedefault_input(:scf)
removedefault_input(:bands)
