#!/bin/sh
# Totally lazy script to run on one of the 72 core workers
# source 3 overloads the memory though so run this one on another worker

JZBIN=$1 # NB: JZBIN denotes source bin number. IZBIN denotes lens bin number. 
echo "Source bin: " $JZBIN

ntheta_bins=4
theta_min=2 #arcmin
theta_max=30 #arcmin
nspins_start=0
nspins_end=500
paramfile=$2     #param_files/params_K1000_BOSS_BlindB_SOMFid.dat



#for JZBIN in 1 2 4 5; do 

for IZBIN in 1 2 3 4 5 
do
    python Shear_ratio_wspin_test.py $ntheta_bins $theta_min $theta_max $JZBIN $nspins_end $IZBIN $nspins_start $paramfile &
    sleep 20s  # give it some time to read the source files before setting off the next run
done    

#done
wait
