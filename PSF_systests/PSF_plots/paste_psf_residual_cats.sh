#!/bin/bash

#patch=$1
allfiles=""

VERSION=V1.0.0
MD=/disk09/KIDS/KIDSCOLLAB_$VERSION/
LFVER=glab_319d
#LFVER=svn_309b

#today=$1
#FIELDLIST=/home/cech/KiDSLenS/THELI_catalogues/ROE_scripts/new_fields/list_$today

patch=$1
FIELDLIST=/home/cech/KiDSLenS/THELI_catalogues/KIDS_conf/$patch.txt

#while read AWfield field dum dum
#do

#    PSFres=$MD/$field/checkplots/PSFRES_XI_$LFVER/${field}*PSFres.cat

#    if [ -d $MD/$field/checkplots/PSFRES_XI_$LFVER ]; then
	allfiles=$allfiles" "$PSFres
#    fi
	
#done < /home/cech/KiDSLenS/THELI_catalogues/KIDS_conf/$patch.txt
#ldacpaste_theli -i $allfiles -o ${patch}_PSF_residuals_${VERSION}_$LFVER.cat

# ldacpaste has a maximum number of files of 1024
# so add in a counter


counter=0
partcounter=0
partfiles=""

#while read field
while read AWfield field dum dum
do

    find $MD/$field/checkplots/*$LFVER -empty -delete
    
    if [ -d $MD/$field/checkplots/PSFRES_XI_$LFVER ]; then

	PSFres=$MD/$field/checkplots/PSFRES_XI_$LFVER/${field}*PSFres.cat
	
	allfiles=$allfiles" "$PSFres
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

#echo $allfiles
if [ $partcounter -gt 0 ]; then
    ldacpaste_theli -i $partfiles -o ${patch}_PSF_residuals_${VERSION}_$LFVER.cat
    rm -f PASTE_TMP_*
else
    ldacpaste_theli -i $allfiles -o ${patch}_PSF_residuals_${VERSION}_$LFVER.cat
fi
