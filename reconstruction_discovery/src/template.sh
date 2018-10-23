#!/usr/bin/env bash
#SBATCH -N 1
#SBATCH -p "clustername"
#SBATCH -A "projectnumber"
#SBATCH -o my.stdout
#SBATCH -t 12:00:00
#SBATCH --mail-user "user@mail"

#Load module in cluster
module load matlab/7.14 

#Export paths
MOSEKPLATFORM=linux64x86
export PATH=$path/to/mosek/7/tools/platform/linux64x86/bin/mosek:$PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/mosek/7/tools/platform/linux64x86/bin
export MOSEKLM_LICENSE_FILE=/path/to/mosek/7/licenses/mosek.lic