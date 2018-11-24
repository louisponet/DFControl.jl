using DFControl, Test
import DFControl: searchdir, QEVariableInfo, QEControlBlockInfo, QEInputInfo, QEDataBlockInfo, qevariable, QEInputInfos

@test qevariable(Exec("projwfc.x"), :calculation).typ == Nothing
@test qevariable(Exec("pw.x"), :calculation).typ == qevariable(:calculation).typ
