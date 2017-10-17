using DFControl, Base.Test

test_job_path = joinpath(@__DIR__,"../assets/inputs/test_job")
df_job = load_job("test_job",test_job_path)
df_job2 = load_job("test_job",joinpath(@__DIR__,test_job_path),new_homedir="blabla")
@test df_job2.home_dir    == "blabla/"
@test length(df_job.calculations) == 6
@test df_job.calculations[3][1] == "~/bin/pw.x"
@test df_job.home_dir     == test_job_path*"/"

mkdir(joinpath(test_job_path,"test_dir/"))
test_dir = joinpath(test_job_path,"test_dir/")
df_job.home_dir = test_dir
save_job(df_job)
df_job2 = load_job("test_job",test_dir)
@test begin
  for (i,(calc_key,calc_tup)) in enumerate(df_job.calculations)
    calc = calc_tup
    calc2 = df_job2.calculations[i][2]
    for (control_key,control_block) in calc.control_blocks
      calc2_control_block = calc2.control_blocks[control_key]
      for (key,value) in control_block
        if typeof(value)<:Array
          continue
        end
        test = isequal(value,calc2_control_block[key])
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

change_data = Dict(:sk1=>3,:sk2=>3.2,:prefix=>"'test'",:noncolin => false, :ecutwfc=> 35,:test => true, :ion_dynamics=>true , :trash=>'d')
change_data2 = Dict(:bleirgh => "'stuff'")
change_job_data!(df_job,change_data)
change_job_data!(df_job,change_data2)
check_keys = Symbol[:sk1,:prefix,:noncolin,:ecutwfc]
@test check_job_data(df_job,check_keys) == Dict(:sk1=>3,:prefix=>"'test'",:noncolin => false, :ecutwfc=> 35)

set_data1 = Dict(:Ze => [Point3D(1.2,3.2,1.2)])
set_data2 = Dict(:control => Dict(:test => true))
set_job_data!(df_job,[1,3],:atoms,set_data1)
set_job_data!(df_job,[1,3],:control_blocks,set_data2)
@test df_job.calculations[3][2].control_blocks[:control][:test]
@test df_job.calculations[1][2].control_blocks[:control][:pseudo_dir] == "'./'"
@test df_job.calculations[1][2].atoms[:Ze] == [Point3D(1.2,3.2,1.2)]

println("")
