using Revise
using DFControl

# change_job_data!(test_job,Dict(:restart => false))
# submit_job(test_job)
test_job = load_job("test_job","assets/inputs/test_job")

test_job.home_dir = "test/test_job"

new_job = load_server_job("test_job","ponet@10.255.9.115","GeTe_2/nonrel/test","test/test_job")
save_job(new_job)
submit_job(new_job)
