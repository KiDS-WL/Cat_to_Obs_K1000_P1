#!/bin/bash -u

# Location of the GGL cats
DD=/disk09/KIDS/K1000_TWO_PT_STATS/GGLCATS/
# The fine binning that we want for the GGL shear ratio test
TOMOINFO="5 0.2 0.3 0.4 0.5 0.6 0.7"
# Extract the number of tomobins from the TOMOINFO
TOMOINFOARR=($TOMOINFO)
NTOMO=${TOMOINFOARR[0]}

for filetype in data random
do


    # Run through the NTOMO tomobins
    for ((i = 1 ; i <= $NTOMO ; i++)); do
        j=$(($i + 1))
        zmin=${TOMOINFOARR[$i]}
        zmax=${TOMOINFOARR[$j]}

	if (( $(echo "$zmax < 0.51" |bc -l) )); then
	    infile=$DD/BOSS_${filetype}_z1.fits
	else
	    infile=$DD/BOSS_${filetype}_z2.fits
	fi

	outfile=$DD/BOSS_${filetype}_${NTOMO}Z_$i.fits
	
	python create_tomocats_GGL.py $zmin $zmax $infile $outfile
    done
done
