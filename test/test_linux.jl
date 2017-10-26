using Revise
using DFControl
using Plots
#TODO visualize the structure
BiAlO3 = load_job("/home/ponet/Documents/PhD/BiAlO3/NSOC")
BiAlO3 = load_server_job("BiAlO3/NSOC","/home/ponet/Documents/PhD/BiAlO3/NSOC")
atoms = get_data(BiAlO3,"bands",:atomic_positions)
atoms[:Te] = [Point3D(0.2f0)]
change_atoms!(BiAlO3,atoms,:pbesol,"paw")