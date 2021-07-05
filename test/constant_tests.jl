using DFControl, Test
import DFControl: searchdir, QEFlagInfo, QEControlBlockInfo, QECalculationInfo,
                  QEDataBlockInfo, qe_flaginfo, QECalculationInfos

@test eltype(qe_flaginfo(Exec(; exec = "projwfc.x"), :calculation)) == Nothing
@test eltype(qe_flaginfo(Exec(; exec = "pw.x"), :calculation)) ==
      eltype(qe_flaginfo(:calculation))
