using Revise
using DFControl

test_wan = read_wannier_input("assets/inputs/wannier/wan.win")
test_wan2 = read_wannier_input("assets/inputs/wannier/wan_test.win")
write_wannier_input("assets/inputs/wannier/wan.wintest",test_wan)
test_job = load_job("test_job","assets/inputs/test_job/",new_homedir="assets/inputs/test_job/test_test/")

test_job.home_dir = "assets/inputs/test_job/test_test"

save_job(test_job)
