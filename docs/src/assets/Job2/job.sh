#!/bin/bash
# Generated by RemoteHPC
# job-name Si

mpirun -np 14  /home/lponet/Software/qe/bin/pw.x  < scf.in > scf.out
mpirun -np 14  /home/lponet/Software/qe/bin/pw.x  < bands.in > bands.out
mpirun -np 14  /home/lponet/Software/qe/bin/pw.x  < nscf.in > nscf.out
mpirun -np 14  /home/lponet/Software/qe/bin/projwfc.x  < projwfc.in > projwfc.out

