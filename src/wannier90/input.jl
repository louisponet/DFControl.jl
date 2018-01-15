"""
    change_k_points!(input::QEInput, k_grid::NTuple{3, Int})

Changes the data in the k point `DataBlock` inside the specified calculation.
"""
function change_k_points!(input::WannierInput, k_grid::NTuple{3, Int})
    change_flags!(input, :mp_grid => [k_grid...])
    k_points = gen_k_grid(k_grid[1], k_grid[2], k_grid[3], :wan)
    change_data!(input, :kpoints, k_points)
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

function change_cell!(input::WannierInput, cell_parameters; option=:ang)
    @assert size(cell_parameters) == (3,3) "Cell parameters has wrong size.\nExpected: (3,3) got ($(size(cell_parameters)[1]), $(size(cell_parameters)[2]))."
    if option == :angstrom
        option = :ang
    end
    change_data!(input, :unit_cell_cart, cell_parameters, option=option)
end

"""
    get_atoms(input::WannierInput)

Returns a list of the atomic positions in Angstrom.
"""
function get_atoms(input::WannierInput)
    out_atoms = OrderedDict{Symbol, Array{<:Point3D, 1}}()
    for block in input.data_blocks
        if contains(string(block.name), "atoms")
            if block.name == :atoms_cart 
                cell = block.option == :bohr ? conversions[:bohr2ang] * eye(3) : eye(3)
            else
                cell_block = get_block(input, :unit_cell_cart)
                alat = cell_block.option == :bohr ? conversions[:bohr2ang] : 1.0
                cell = alat * cell_block.data 
            end
            for (atom, positions) in block.data
                t_pos = Point3D[]
                for pos in positions
                    push!(t_pos, cell' * pos) 
                end
                out_atoms[atom] = t_pos
            end
        end
    end
    return out_atoms
end

"""
    change_atoms!(input::WannierInput, atoms::OrderedDict{Symbol,<:Array{<:Point3D,1}}; option=:angstrom)

Changes the atoms in the input to the specified atoms.
"""
function change_atoms!(input::WannierInput, atoms::OrderedDict{Symbol,<:Array{<:Point3D,1}}; option=:ang, kwargs...)
    if option == :angstrom
        option = :ang
    end
    for block in input.data_blocks
        if  block.name == :atoms_frac
            block.name = :atoms_cart
            change_data!(input, :atoms_cart, atoms, option=option)
        elseif block.name == :atoms_cart
            change_data!(input, :atoms_cart, atoms, option=option)
        end
    end
end


