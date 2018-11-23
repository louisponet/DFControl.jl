using DFControl, Test
import DFControl: searchdir, QEVariableInfo, QEControlBlockInfo, QEInputInfo, QEDataBlockInfo, qevariable, QEInputInfos

@test qevariable(QEInputInfos[2], :calculation).typ == Nothing
@test qevariable(QEInputInfos[3], :calculation).typ == qevariable(:calculation).typ
