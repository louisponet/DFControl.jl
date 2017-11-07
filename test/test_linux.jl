using Revise
using DFControl
using Plots
#TODO visualize the structure
BiAlO3 = load_job("/home/ponet/Documents/PhD/BiAlO3/NSOC")
BiAlO3 = load_server_job("BiAlO3/NSOC","/home/ponet/Documents/PhD/BiAlO3/NSOC")
pull_outputs(BiAlO3, extras = ["*.xsf"])
Pkg.test("DFControl")

