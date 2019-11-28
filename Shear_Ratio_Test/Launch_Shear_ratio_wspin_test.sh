#!/bin/bash

#SBATCH -n 2
#SBATCH -p all
#SBATCH -t 7-00:00
#SBATCH --job-name=Shear_Spin_Ratio
#SBATCH --requeue
#SBATCH --mail-type=ALL
#SBATCH --constraint="3TBdatadisk"
#SBATCH --mem=300000

JZBIN=$1
###### NB: Normally -mem=150000

###################### mem is in MB

./run_shear_ratio_w_spin.sh $JZBIN


