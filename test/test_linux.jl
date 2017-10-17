using Revise
using DFControl

test_job = load_qe_job("Test_job","assets/inputs/qe/")
test_job.home_dir="assets/inputs/qe/tesjob/"
save_job(test_job)
test_wan = read_wannier_input("assets/inputs/wannier/wan.win")
