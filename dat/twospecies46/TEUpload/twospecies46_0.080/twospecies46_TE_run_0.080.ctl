#!/bin/bash
# parallel job using 48 processors.  Runs for 6.0 hours (max).
#PBS -l nodes=6:ppn=8,walltime=6:00:00
#PBS -m abe
#PBS -M chaneyl@princeton.edu
module load openmpi/gcc/1.4.5/64
module load hdf5/gcc/openmpi-1.4.5/1.8.8
cd /home/chaneyl/twospecies46/TE/twospecies46_0.080
mkdir out
mpiexec -n 48 /home/florescu/Exec/MPI_MFT_v1/bin/mft-mpi \
mesh-size=7 \
epsAir=1 epsSi=11.56 \
CylinderHorizontalRadius=0.06 \
num-bands=650 \
k-interp=3 res=24 res_x=12 res_y=12 \
twospecies46_TE_param_0.080.ctl > ./out/twospecies46_TE_out_0.080.OUT ; 
