using Revise
using DFControl
using Plots
#TODO visualize the structure
server = "ponet@10.255.9.115"
job = load_job("/home/ponet/Documents/PhD/BiTeI/NSOC/",server=server,server_dir="BiTeI/NSOC");
change_atoms!(job,get_data(read_qe_input("assets/inputs/qe/scf.in"),:atomic_positions))
print_flow(job)
print_info(job)
set_should_run!(job,[(2,true),(3,true),(4,true)])
outputs = pull_outputs(job,output_fuzzy=["*_hr.dat","*_r.dat"])
bands = read_qe_bands_file(job.local_dir*"2BiTeI_bands.out")
get_flag(job,:dis_win_min)
plot(bands,ylims=[0,15])
submit_job(job)


print_run_command(job,"nscf")