#!/bin/bash
# script to combine all PSF residual files into a single mosaic file
# uses ldacpaste

allfiles=""

VERSION=V1.0.0
MD=/disk09/KIDS/KIDSCOLLAB_$VERSION/
LFVER=glab_321
#LFVER=svn_309b

patch=$1
FIELDLIST=/home/cech/KiDSLenS/THELI_catalogues/KIDS_conf/$patch.txt

# ldacpaste has a maximum number of files of 1024
# so add in a counter

counter=0
partcounter=0
partfiles=""

while read field
#while read AWfield field dum dum
do
    # check that the residual files directory exists  
    if [ -d $MD/$field/checkplots/PSFRES_XI_$LFVER ]; then

	# check that files are contained within it

	nfiles=`ls -l $MD/$field/checkplots/PSFRES_XI_$LFVER/${field}*PSFres.cat | wc -l`

	if [ "$nfiles" -eq "0" ]; then
	    echo "$nfiles no PSF checkplot files for field $field"
	else
	
	    PSFres=$MD/$field/checkplots/PSFRES_XI_$LFVER/${field}*PSFres.cat
	    allfiles=$allfiles" "$PSFres

	fi
	    
    else
	echo "Missing $MD/$field/checkplots/PSFRES_XI_$LFVER"
	exit
    fi

    let counter=counter+1

    if [ $counter -gt 100 ]; then
	    counter=0
	    let partcounter=partcounter+1
	
	    ldacpaste_theli -i $allfiles -o PASTE_TMP_$partcounter

        allfiles=""

	    partfiles=$partfiles" "PASTE_TMP_$partcounter
    fi
	
done < $FIELDLIST

if [ $partcounter -gt 0 ]; then
    ldacpaste_theli -i $partfiles -o ${patch}_PSF_residuals_${VERSION}_$LFVER.cat
    rm -f PASTE_TMP_*
else
    ldacpaste_theli -i $allfiles -o ${patch}_PSF_residuals_${VERSION}_$LFVER.cat
fi
