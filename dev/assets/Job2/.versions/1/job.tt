#!/bin/bash
#SBATCH -J Si 

#mpirun -np 4 /opt/qe/bin/pw.x -nk 4 < scf.in > scf.out
mpirun -np 4 /opt/qe/bin/pw.x -nk 4 < bands.in > bands.out
mpirun -np 4 /opt/qe/bin/pw.x -nk 4 < nscf.in > nscf.out
mpirun -np 4 /opt/qe/bin/projwfc.x < projwfc.in > projwfc.out
