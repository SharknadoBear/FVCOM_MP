#!/bin/bash

module purge

# Updated module loading to match Lmod automatic replacements
# Removed gcc/12.1.0 since Lmod automatically replaces it with intel/2023.2.0
module load intel/2023.2.0
module load intel-oneapi-mpi/2021.11.0-intel
module load netcdf-c/4.9.2-cray-mpich-intel
module load netcdf-fortran/4.6.1-intel

# Display loaded modules for verification
echo "Loaded modules:"
module list