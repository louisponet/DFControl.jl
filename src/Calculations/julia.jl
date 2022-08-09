struct Julia <: Package end

const JULIA_EXECS = "julia"

flagtype(::Calculation{Julia}, k) = k in (:source, :script) ? String : Nothing
