# using GLVisualize, GLWindow, GLAbstraction,Colors,GeometryTypes
using  GtkReactive, Gtk.ShortNames, Colors,Gtk
# const dash_window = glscreen("DFDashboard",resolution=(1600,1200))
const dash_window = Window("DFDashboard",400,800)
grid = Grid()
bx = Box(:v)
scwindow1 = ScrolledWindow()
adj = getproperty(scwindow1,:vadjustment, GtkAdjustment)


# scwindow2 = ScrolledWindow()
setproperty!(scwindow1,:height_request,800)
setproperty!(scwindow1,:width_request,400)
setproperty!(scwindow1,:vexpand,true)
# setproperty!(scwindow2,:height_request,800)
text_a = Box(:h)
# text_b = Box(:h)
push!(scwindow1,text_a)
# push!(scwindow2,text_b)
grid[1,1] = scwindow1
# grid[2,1] = scwindow2
setproperty!(grid,:column_homogeneous,true)
setproperty!(grid,:column_spacing,15)
push!(dash_window,grid)
# push!(dash_window,scwindow2)
acc(x,y) = x*"\n"*y
tb = textarea(signal=foldp(acc,"",print_s))
function scroll_cb(widgetptr::Ptr, rectptr::Ptr, user_data)
  setproperty!(adj,:value, getproperty(adj,:upper,AbstractFloat) - getproperty(adj,:page_size,AbstractFloat))
  nothing
end
signal_connect(scroll_cb, tb.widget, "size-allocate", Void, (Ptr{Gtk.GdkRectangle},), false, ())
# tb2 = textarea(signal=foldp(acc,"",print_s))
setproperty!(tb,:editable,false)
# setproperty!(tb2,:editable,false)
setproperty!(tb,:can_focus,false)
# setproperty!(tb2,:can_focus,false)
push!(text_a,tb)
# push!(text_a,tb2)
showall(dash_window)

const colors = Dict{Symbol,RGBA}()

colors[:text_background] = RGBA(24/255,32/255,38/255,1.f0)
colors[:text_color]      = RGBA(1.0f0,1.0f0,1.0f0,1.0f0)

# function init_dash(w::Screen)
#   text, visa = x_partition(w.area, 50)
#   texta,textb = y_partition(text,50)

#   text_screen   = Screen(w, name=:text,area=text)
#   text_screen_a = Screen(text_screen, name=:texta, area = texta, color=colors[:text_background])
#   text_screen_b = Screen(text_screen, name=:textb, area = textb, color=colors[:text_background]) 

#   # text_b = visualize(preserve(foldp(*,"s",preserve(droprepeats(print_s)))),color=colors[:text_color],scale=Vec2f0(50),stroke_width=1f0,scale_primitive=false)
#   # interval = preserve(every(5.))
#   acc(x,y) = x*"\n"*y
#   text_accum = foldp(acc,"t",print_s)
#   preserve(map(println,map(length,text_accum)))
#   # text_to_screen = sampleon(interval,text_accum)
#   text_b = visualize(text_accum,color=colors[:text_color])
#   _view(text_b,text_screen_b,camera=:orthographic_pixel)
#   center!(text_screen_b,:orthographic_pixel)

#   @async begin
#     renderloop(w)
#   end
# end

# init_dash(dash_window)

# function reload_dash()
#   dash_window = glscreen("DFDashboard",resolution=(1600,1200))
#   init_dash(dash_window)
# end
