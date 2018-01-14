using Revise
using DFControl
using Plots

#TODO visualize the structure
reload_dash()
BiAlO3 = load_job("/home/ponet/Documents/PhD/BiAlO3/NSOC")
print_flags(BiAlO3)
print_flags(BiAlO3)
BiAlO3 = load_server_job("BiAlO3/NSOC","/home/ponet/Documents/PhD/BiAlO3/NSOC")
pull_outputs(BiAlO3, extras = ["*.xsf"])
Pkg.test("DFControl")

using GtkReactive, Gtk.ShortNames
const dash_window = Window("DFDashboard",1800,1200) |> (bx = Box(:v));
push!(dash_window,bx)
tb = textbox(String)
push!(bx,tb)
showall(dash_window)
push!(tb,"test")

win = Window("Testing") |> (bx = Box(:v));
sl = slider(1:11)
push!(bx, sl);
showall(win);
using Gtk
scwindow=ScrolledWindow()


