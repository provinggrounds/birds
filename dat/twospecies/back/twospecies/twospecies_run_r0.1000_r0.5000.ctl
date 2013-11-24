#!/bin/bash
# parallel job using 48 processors.  Runs for 2.5 hours (max).
#PBS -l nodes=6:ppn=8,walltime=2:30:00
#PBS -m abe
#PBS -M chaneyl@princeton.edu
module load openmpi/gcc/1.4.5/64
module load hdf5/gcc/openmpi-1.4.5/1.8.8
cd /home/chaneyl/twospecies
mkdir out
mpiexec -n 48 /home/florescu/Exec/MPI_MFT/bin/mft-mpi \
mesh-size=7 \
epsAir=1 epsSi=11.56 \
CylinderVerticalRadius=0.189  \
num-bands=350 \
k-interp=3 res=24 res_x=12 res_y=12 \
twospecies_param_r0.1000_r0.5000.ctl > ./out/twospecies_out_r0.1000_r0.5000.OUT ; 
