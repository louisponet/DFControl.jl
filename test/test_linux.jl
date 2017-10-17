using Revise
using DFControl

test_job = load_server_job("test_job","ponet@10.255.9.115","GeTe_2/nonrel/test","test/test_job")
submit_job(test_job)
