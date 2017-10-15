using DFControl, Base.Test
df_job = load_qe_job("test_job",joinpath(@__DIR__,"../assets/inputs/qe"))

@test length(df_job.flow) == 4
@test df_job.flow[3]      == ("~/bin/pw.x","bands.in")
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

change_data = Dict(:prefix=>"'test'",:noncolin => false, :ecutwfc=> 35)
change_job_data!(df_job,change_data)
check_keys = keys(change_data)
@test check_job_data(df_job,check_keys) == change_data
