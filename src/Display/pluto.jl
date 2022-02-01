# using ANSIColoredPrinters

function color_print(io, v)
    tio = IOBuffer()
    tio_ =  IOContext(tio, :color=>true)
    df_show(tio_, v)
    html_str = sprint(io2->display(io2, MIME"text/html"(),
                      HTMLPrinter(tio_, root_class="documenter-example-output")))
    print(io,"$html_str")
end

import Base: display 

display(io::IO, ::MIME"text/html", v::InputData) = color_print(io, v)

display(io::IO, ::MIME"text/html", v::Band) = color_print(io, v)

display(io::IO, ::MIME"text/html", v::Job) = color_print(io, v)
display(io::IO, ::MIME"text/html", v::Calculation) = color_print(io, v)
display(io::IO, ::MIME"text/html", v::Calculations.QEFlagInfo) = color_print(io, v)
display(io::IO, ::MIME"text/html", v::Calculations.ElkFlagInfo) = color_print(io, v)
display(io::IO, ::MIME"text/html", v::Calculations.ElkControlBlockInfo) = color_print(io, v)
display(io::IO, ::MIME"text/html", v::Element) = color_print(io, v)
display(io::IO, ::MIME"text/html", v::Structure) = color_print(io, v)
display(io::IO, ::MIME"text/html", v::Atom) = color_print(io, v)
display(io::IO, ::MIME"text/html", v::Projection) = color_print(io, v)

display(io::IO,  v::InputData) = df_show(io, v)

display(io::IO,  v::Band) = df_show(io, v)

display(io::IO,  v::Job) = df_show(io, v)
display(io::IO,  v::Calculation) = df_show(io, v)
display(io::IO,  v::Calculations.QEFlagInfo) = df_show(io, v)
display(io::IO,  v::Calculations.ElkFlagInfo) = df_show(io, v)
display(io::IO,  v::Calculations.ElkControlBlockInfo) = df_show(io, v)
display(io::IO,  v::Element) = df_show(io, v)
display(io::IO,  v::Structure) = df_show(io, v)
display(io::IO,  v::Atom) = df_show(io, v)
display(io::IO,  v::Projection) = df_show(io, v)
