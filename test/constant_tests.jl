using DFControl, Test
import DFControl: searchdir, QEFlagInfo, QEControlBlockInfo, QEInputInfo, QEDataBlockInfo, qe_flaginfo, QEInputInfos

@test qe_flaginfo(Exec("projwfc.x"), :calculation).typ == Nothing
@test qe_flaginfo(Exec("pw.x"), :calculation).typ == qe_flaginfo(:calculation).typ
