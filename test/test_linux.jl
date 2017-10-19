using Revise
using DFControl
new_job = load_job("test/test_job",server = "ponet@10.255.9.115", server_dir = "GeTe_2/nonrel/test")
remove_job_control_flag!(new_job,[6,8], :restart)
set_should_run!(new_job,[false,false,false,false,false,true,true,true])
submit_job(new_job)
#TODO make change and set job data work again
