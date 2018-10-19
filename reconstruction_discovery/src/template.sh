#!/usr/bin/env bash
#SBATCH -N 1
#SBATCH -p sysbio
#SBATCH -A C3SE999-12-7
#SBATCH -o my.stdout
#SBATCH -t 12:00:00
#SBATCH --mail-user gatto@chalmers.se

module load matlab/7.14

MOSEKPLATFORM=linux64x86
export PATH=$/c3se/users/gatto/Glenn//mosek/7/tools/platform/linux64x86/bin/mosek:$PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/c3se/users/gatto/Glenn/mosek/7/tools/platform/linux64x86/bin
export MOSEKLM_LICENSE_FILE=/c3se/users/gatto/Glenn/mosek/7/licenses/mosek.lic
