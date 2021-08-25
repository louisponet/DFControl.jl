#printing that is not needed in Atom
import Base: show

show(io::IO, block::InputData) = df_show(io, block)

show(io::IO, band::Band) = df_show(io, band)

show(io::IO, job::Job) = df_show(io, job)
show(io::IO, c::Calculation) = df_show(io, c)
show(io::IO, flag_info::Calculations.QEFlagInfo) = df_show(io, flag_info)
show(io::IO, flag_info::Calculations.ElkFlagInfo) = df_show(io, flag_info)
show(io::IO, info::Calculations.ElkControlBlockInfo) = df_show(io, info)
show(io::IO, el::Element) = df_show(io, el)
show(io::IO, str::Structure) = df_show(io, str)
show(io::IO, at::Atom) = df_show(io, at)
show(io::IO, proj::Projection) = df_show(io, proj)
