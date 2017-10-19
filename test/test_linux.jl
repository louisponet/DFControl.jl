using Revise
using DFControl
using Plots
#TODO visualize the structure
server = "ponet@10.255.9.115"
job = load_server_job(server,"BiTeI/NSOC/","/home/ponet/Documents/PhD/BiTeI/NSOC")

job.calculations[5].data_blocks[end].data = [k[1:3] for k in job.calculations[4].data_blocks[end].data]
job.calculations[4].run = true
set_should_run!(job,[false,true,false,true,true,true,true])
remove_job_control_flags!(job,[:occupations])
submit_job(job)



i
