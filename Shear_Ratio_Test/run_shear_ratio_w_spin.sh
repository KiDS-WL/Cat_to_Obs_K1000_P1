#!/bin/sh
# Totally lazy script to run on one of the 72 core workers
# source 3 overloads the memory though so run this one on another worker

ntheta_bins=4
theta_min=2 #arcmin
theta_max=30 #arcmin
nspins=1000

#for JZBIN in 3 #1 2 3 4 5 
JZBIN=3
for IZBIN in 4 5 #3 #4 5 
do  
    python Shear_ratio_wspin_test.py $ntheta_bins $theta_min $theta_max $JZBIN $nspins $IZBIN 600 &
    sleep 20s  # give it some time to read the source files before setting off the next run
done    

wait