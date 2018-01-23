"""
    change_kpoints!(input::QEInput, k_grid::NTuple{3, Int})

Changes the data in the k point `DataBlock` inside the specified calculation.
"""
function change_kpoints!(input::WannierInput, k_grid::NTuple{3, Int}; print=true)
    change_flags!(input, :mp_grid => [k_grid...])
    k_points = gen_k_grid(k_grid[1], k_grid[2], k_grid[3], :wan)
    change_data!(input, :kpoints, k_points, print=print)
end

"Returns the cell parameters in Angstrom."
function get_cell(input::WannierInput)
    block = get_block(input, :unit_cell_cart)
    if block != nothing
        if block.option == :bohr
            alat = conversions[:bohr2ang]
        else
            alat = 1.0
        end
        return block.data * alat
    end
end

function change_cell!(input::WannierInput, cell_parameters; option=:ang, print=true)
    @assert size(cell_parameters) == (3,3) "Cell parameters has wrong size.\nExpected: (3,3) got ($(size(cell_parameters)[1]), $(size(cell_parameters)[2]))."
    if option == :angstrom
        option = :ang
    end
    change_data!(input, :unit_cell_cart, cell_parameters, option=option, print=print)
end