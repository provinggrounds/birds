#!/bin/bash
# parallel job using 8 nodes and 16 processors.  Runs for 4 hours (max).
#PBS -l procs=48,walltime=2:30:00
module load openmpi/gcc/1.4.5/64
module load hdf5/gcc/openmpi-1.4.5/1.8.8
cd /home/chaneyl/twospecies/twospecies_radtest
mpiexec -n 48 /home/florescu/Exec/MPI_MFT/bin/mft-mpi \
mesh-size=7 \
epsAir=1 epsSi=11.56 \
CylinderVerticalRadius=0.189  \
num-bands=600 \
k-interp=3 res=24 res_x=12 res_y=12 \
spT_irad0.022_rat15.ctl > spT_irad0.022_rat15.out ; 
