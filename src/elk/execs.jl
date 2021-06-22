const ELK_EXECS = ["elk", "elk-omp"]

is_elk_exec(exec::Exec) = exec.exec ∈ ELK_EXECS
