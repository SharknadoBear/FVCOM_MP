#!/bin/bash

#SBATCH --account=hindcastra
#SBATCH --time=02-00:00:00
#SBATCH --job-name=WFPrun
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=104
#SBATCH -o OUTPUT.%j
#SBATCH -e ERROR.%j
#SBATCH  --partition=standard
##SBATCH --mail-type=begin
##SBATCH --mail-type=end
#SBATCH --mail-user=yicheng.huang@pnnl.gov


module purge

##gcc
#module load gcc/13.1.0
#module load hdf5/1.14.1-2-openmpi-gcc
#module load openmpi/4.1.1-gcc
#module load netcdf-c/4.9.2-openmpi-gcc
#module load netcdf-fortran/4.6.0-gcc

##intel
#module load intel/2023.2.0
#module load mpich/4.1-intel
#module load hdf5/1.14.1-2-intel-oneapi-mpi-intel
#module load netcdf-c/4.9.2-intel-oneapi-mpi-intel
#module load netcdf-fortran/4.6.0-intel

##intel-cray
module load craype-x86-spr
module load gcc/12.1.0
module load PrgEnv-intel
module swap cray-mpich cray-mpich-abi
module load netcdf-c/4.9.2-cray-mpich-intel
module load netcdf-fortran/4.6.1-intel
module load hdf5/1.14.1-2-intel-oneapi-mpi-intel

export NETCDF=/nopt/nrel/apps/cpu_stack/libraries-craympich/06-24/linux-rhel8-sapphirerapids/oneapi-2023.2.0/netcdf-fortran-4.6.1-dzwy5uqubyodxnvwyj4pimnkm63tx6qy
export WRFIO_NCD_LARGE_FILE_SUPPORT=1
export NETCDF_classic=1
export OMP_NUM_THREADS=1
export OMP_PROC_BIND=spread
export KMP_AFFINITY=balanced

export FI_CXI_RX_MATCH_MODE=hybrid


ulimit -s unlimited

echo "Start:" `date`

#mpirun -np 528 ./coawstM coupling_WFIP3.in
srun ./fvcom_TLMF --casename=waterPACT_c

