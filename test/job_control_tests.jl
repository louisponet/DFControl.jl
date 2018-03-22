using DFControl, Base.Test

test_job_path = joinpath(@__DIR__, "test_job")
df_job = load_job(test_job_path);
df_job2 = load_job(joinpath(@__DIR__, test_job_path), new_local_dir="blabla");
@test df_job2.local_dir    == "blabla/"
@test length(df_job.calculations) == 7
@test input(df_job, "nscf").runcommand.exec == runcommand(df_job, "nscf").exec
@test df_job.local_dir     == test_job_path*"/"

try mkdir(joinpath(test_job_path,"test_dir/")) end
test_dir = joinpath(test_job_path,"test_dir/")
df_job.local_dir = test_dir
save(df_job)
df_job2 = load_job(test_dir)
@test begin
  for (i,calc) in enumerate(df_job.calculations)
    calc2 = df_job2.calculations[i]
    for (j,block) in enumerate(calc.data)
      block2 = calc2.data[j]
      for name in fieldnames(block)
        field1 = getfield(block,name)
        field2 = getfield(block2,name)
        test = isequal(field1,field2)
        if !test
          return false
        end
      end
    end
  end
  return true
end

files_to_remove = DFControl.search_dir(test_dir,".")
for file in files_to_remove
  rm(test_dir * file)
end
rm(test_dir)

# setdata = Dict(:sk1=>3,:sk2=>3.2,:prefix=>"'test'",:noncolin => false, :ecutwfc=> 35,:test => true, :ion_dynamics=>true , :trash=>'d')
# setdata2 = Dict(:bleirgh => "'stuff'")

scf_input     = read_qe_input(joinpath(test_job_path, "scf.in"))[1]
t_l = length(df_job.calculations)
add!(df_job, scf_input)
@test length(df_job.calculations)==t_l+1
setflags!(df_job, :dis_win_min => 9.2)
@test flag(df_job, "wan.win", :dis_win_min) == 9.2
setflags!(df_job, "nscf", :calculation => "'scf'")
@test flag(df_job,"nscf",:calculation)== "'scf'"

setflags!(df_job, :Hubbard_J => 1.0)
@test flag(df_job, "nscf", :Hubbard_J) == 1.0

setflags!(df_job, :dis_win_min => 24.0)
@test flag(df_job, "wan.win", :dis_win_min) == 24.0
rmflags!(df_job,"wan.win", :dis_win_min, :dis_win_max)

setflow!(df_job, [false for i=1:length(df_job.calculations)])
@test df_job.calculations[1].run == false
setflow!(df_job, "nscf" => true, "bands" => true)
@test df_job.calculations[3].run
setflow!(df_job, "pw2wan" => true)
@test df_job.calculations[end-1].run

runcommand!(df_job, "nscf", Exec("test"))
@test runcommand(df_job, "nscf").exec == Exec("test").exec

print_run_command(df_job, "nscf")
print_flow(df_job)
@test print_block(df_job,:atomic_positions) == print_block(df_job,"nscf",:atomic_positions)
@test print_blocks(df_job) == print_data(df_job)
@test print_blocks(df_job,"bands") == print_data(df_job,"bands")
print_info(df_job)
print_flags(df_job)
@test print_flags(df_job,"nscf") == print_flags(df_job,["nscf"])
@test print_flags(df_job,[:dis_win_min]) == print_flag(df_job,:dis_win_min)
@test inputs(df_job,["nscf"]) == inputs(df_job,"nscf")







# @test check_job_data(df_job,check_keys) == Dict(:sk1=>3,:prefix=>"'test'",:noncolin => false, :ecutwfc=> 35)

# setdata1 = Dict(:Ze => [Point3(1.2,3.2,1.2)])
# setdata2 = Dict(:control => Dict(:test => true))
# setjob_data!(df_job,[1,3],:atoms,setdata1)
# setjob_data!(df_job,[1,3],:control,setdata2)
# @test df_job.calculations[3][2].control[:control][:test]
# @test df_job.calculations[1][2].control[:control][:pseudo_dir] == "'./'"
# @test df_job.calculations[1][2].atoms[:Ze] == [Point3(1.2,3.2,1.2)]

# dfprintln("")
