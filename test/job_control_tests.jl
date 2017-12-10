using DFControl, Base.Test

test_job_path = joinpath(@__DIR__,"../assets/inputs/test_job")
df_job = load_job(test_job_path)
display(df_job)
df_job2 = load_job(joinpath(@__DIR__,test_job_path),new_local_dir="blabla")
@test df_job2.local_dir    == "blabla/"
@test length(df_job.calculations) == 7
@test get_input(df_job,"nscf").run_command == get_run_command(df_job,"nscf")
@test df_job.local_dir     == test_job_path*"/"

try mkdir(joinpath(test_job_path,"test_dir/")) end
test_dir = joinpath(test_job_path,"test_dir/")
df_job.local_dir = test_dir
save_job(df_job)
df_job2 = load_job(test_dir)
@test begin
  for (i,calc) in enumerate(df_job.calculations)
    calc2 = df_job2.calculations[i]
    for (j,block) in enumerate(calc.data_blocks)
      block2 = calc2.data_blocks[j]
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

# change_data = Dict(:sk1=>3,:sk2=>3.2,:prefix=>"'test'",:noncolin => false, :ecutwfc=> 35,:test => true, :ion_dynamics=>true , :trash=>'d')
# change_data2 = Dict(:bleirgh => "'stuff'")
data = get_data(df_job,"nscf",:atomic_positions)
data[:Te] = data[:Te].+Point3D(0.01f0)
change_data!(df_job,"nscf",:atomic_positions,data)
@test get_block(df_job,"nscf",:atomic_positions).data == data

scf_input     = read_qe_input(joinpath(@__DIR__,"../assets/inputs/qe/scf.in"),Float64)
t_l = length(df_job.calculations)
add_calculation!(df_job, scf_input)
@test length(df_job.calculations)==t_l+1
change_flags!(df_job,Dict(:dis_win_min => 9.2f0))
@test get_flag(df_job,"wan.win",:dis_win_min)==9.2f0
change_flags!(df_job,["nscf"],Dict(:calculation => "'scf'"))
@test get_flag(df_job,"nscf",:calculation)== "'scf'"

add_flags!(df_job,:control,Dict(:test => "test"))
@test get_flag(df_job,"nscf",:test) == "test"

add_flags!(df_job,Dict(:test=>"test"))
@test get_flag(df_job,"wan.win",:test)=="test"
remove_flags!(df_job,"wan.win",[:test,:dis_win_max])
remove_flags!(df_job,:test)
@test get_flag(df_job,:test)==nothing
@test get_flag(df_job,:dis_win_max)==nothing

set_flow!(df_job,[false for i=1:length(df_job.calculations)])
@test df_job.calculations[1].run == false
change_flow!(df_job,[(1,true)])
@test df_job.calculations[1].run
change_flow!(df_job,["nscf","bands"],true)
@test df_job.calculations[2].run
change_flow!(df_job,[("pw2wan",true)])
@test df_job.calculations[end-2].run

change_run_command!(df_job,"nscf","test")
@test get_run_command(df_job,"nscf") == "test"

print_run_command(df_job,"nscf")
print_flow(df_job)
@test print_block(df_job,:atomic_positions) == print_block(df_job,"nscf",:atomic_positions)
@test print_blocks(df_job) == print_data(df_job)
@test print_blocks(df_job,"bands") == print_data(df_job,"bands")
print_info(df_job)
print_flags(df_job)
@test print_flags(df_job,"nscf") == print_flags(df_job,["nscf"])
@test print_flags(df_job,[:dis_win_min]) == print_flag(df_job,:dis_win_min) 
@test get_inputs(df_job,["nscf"]) == get_inputs(df_job,"nscf")







# @test check_job_data(df_job,check_keys) == Dict(:sk1=>3,:prefix=>"'test'",:noncolin => false, :ecutwfc=> 35)

# set_data1 = Dict(:Ze => [Point3D(1.2,3.2,1.2)])
# set_data2 = Dict(:control => Dict(:test => true))
# set_job_data!(df_job,[1,3],:atoms,set_data1)
# set_job_data!(df_job,[1,3],:control_blocks,set_data2)
# @test df_job.calculations[3][2].control_blocks[:control][:test]
# @test df_job.calculations[1][2].control_blocks[:control][:pseudo_dir] == "'./'"
# @test df_job.calculations[1][2].atoms[:Ze] == [Point3D(1.2,3.2,1.2)]

# dfprintln("")
