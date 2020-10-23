#!/bin/bash

#SBATCH -n 1
#SBATCH -p all
#SBATCH -t 7-00:00
#SBATCH --job-name=Shear_Spin_Ratio
#SBATCH --requeue
#SBATCH --mail-type=ALL
#SBATCH --constraint="3TBdatadisk"
#SBATCH --mem=150000
JZBIN=$1
paramfile=param_files/params_K1000_GAMA_BlindC_SOMFid.dat

./run_shear_ratio_w_spin.sh $JZBIN $paramfile



