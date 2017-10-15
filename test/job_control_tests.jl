using DFControl, Base.Test

df_job = load_qe_job("test_job",joinpath(@__DIR__,"../assets/inputs/qe"))

@test length(df_job.flow) == 4
@test df_job.flow[3]      == ("~/bin/pw.x","bands.in")
@test df_job.home_dir     == joinpath(@__DIR__,"../assets/inputs/qe/")