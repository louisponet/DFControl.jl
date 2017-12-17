using DFControl
using PyCall

try
    parse("5")
catch
    "5"
end


@pyimport abipy.lessons.lesson_dos_bands as lesson_dos_bands
lesson = lesson_dos.Lesson()
flow = lesson[:make_flow()


#maybe use .attached inputs 