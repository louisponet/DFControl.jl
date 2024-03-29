#!/bin/bash
#SBATCH --job-name=Ni
export OMP_NUM_THREADS=1
mpirun -np 4 ~/Software/qe/bin/pw.x -nk 4 < vcrelax.in > vcrelax.out
mpirun -np 4 ~/Software/qe/bin/pw.x -nk 4 < scf.in > scf.out
mpirun -np 4 ~/Software/qe/bin/pw.x -nk 4 < bands.in > bands.out
mpirun -np 4 ~/Software/qe/bin/pw.x -nk 4 < nscf.in > nscf.out
mpirun -np 4 ~/Software/qe/bin/projwfc.x < projwfc.in > projwfc.out
~/Software/wannier90/wannier90.x -pp wanup.win > wanup.wout
mpirun -np 4 ~/Software/qe/bin/pw2wannier90.x < pw2wan_wanup.in > pw2wan_wanup.out
~/Software/wannier90/wannier90.x wanup.win > wanup.wout
~/Software/wannier90/wannier90.x -pp wandn.win > wandn.wout
mpirun -np 4 ~/Software/qe/bin/pw2wannier90.x < pw2wan_wandn.in > pw2wan_wandn.out
~/Software/wannier90/wannier90.x wandn.win > wandn.wout
