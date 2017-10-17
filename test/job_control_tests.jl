using DFControl, Base.Test


df_job = load_qe_job("test_job",joinpath(@__DIR__,"../assets/inputs/qe"))
df_job2 = load_qe_job("test_job",joinpath(@__DIR__,"../assets/inputs/qe"),new_homedir="blabla")
@test df_job2.home_dir    == "blabla/"
@test length(df_job.flow) == 4
@test df_job.flow[3]      == ("~/bin/pw.x","bands")
@test df_job.home_dir     == joinpath(@__DIR__,"../assets/inputs/qe/")

mkdir(joinpath(@__DIR__,"../assets/inputs/qe/test_dir/"))
test_dir = joinpath(@__DIR__,"../assets/inputs/qe/test_dir/") 
df_job.home_dir = test_dir 
save_job(df_job)
df_job2 = load_qe_job("test_job",test_dir)
@test begin 
  for (calc_key,calc) in df_job.calculations
    calc2 = df_job2.calculations[calc_key]
    for (control_key,control_block) in calc.control_blocks
      calc2_control_block = calc2.control_blocks[control_key]
      for (key,value) in control_block 
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

change_data = Dict(:sk1=>3,:sk2=>3.2,:prefix=>"'test'",:noncolin => false, :ecutwfc=> 35,:test => true, :ion_dynamics=>true , :kaka=>'d')
change_data2 = Dict(:bleirgh => "'stuff'")
change_job_data!(df_job,change_data)
change_job_data!(df_job,change_data2)
check_keys = Symbol[:sk1,:prefix,:noncolin,:ecutwfc] 
@test check_job_data(df_job,check_keys) == Dict(:sk1=>3,:prefix=>"'test'",:noncolin => false, :ecutwfc=> 35)

set_data1 = Dict(:Ze => [Point3D(1.2,3.2,1.2)])
set_data2 = Dict(:control => Dict(:test => true))
set_job_data!(df_job,["bands","scf"],:atoms,set_data1)
set_job_data!(df_job,["bands","scf"],:control_blocks,set_data2)
@test df_job.calculations["bands"].control_blocks[:control][:test]
@test df_job.calculations["scf"].control_blocks[:control][:pseudo_dir] == "'./'"
@test df_job.calculations["scf"].atoms[:Ze] == [Point3D(1.2,3.2,1.2)]

println("")