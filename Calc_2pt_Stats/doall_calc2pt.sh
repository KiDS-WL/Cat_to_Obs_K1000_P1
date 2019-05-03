#!/bin/bash -u

# ----------------------------------------------------------------
# File Name:           doall_calc2pt.sh
# Author:              Marika Asgari (ma@roe.ac.uk)
#                      Catherine Heymans (heymans@roe.ac.uk)
# Description:         Master script to do all KiDS tasks related to converting a catalogue
#                      to the cosmic shear and GGL observables that we want to measure
# ----------------------------------------------------------------

# Script history information:
# 2nd May 2019:  Started based on do_all script for pipeline processing - thanks to Thomas and co

##
## function definitions:
##
function printUsage
{
  echo "SCRIPT NAME:"
  echo "    doall_calc2pt.sh"
  echo ""
  echo "    Master 'doall' catalogue processing script to create observables"
  echo ""
  echo "ARGUMENTS:"
  echo "    must be passed with the following switches. Only the mode switch"
  echo "    is necessary. If other arguments are not supplied, the script"
  echo "    will work from default settings. Possible arguments are:"
  echo "       -d /path/to/catalogues"
  echo "       -o /path/to/results"
  echo "       -p patch name N or S"
  echo "       -m list of modes"
  echo "       -v lensfit version"
  echo "       -n ntomo number of tomographic source bins, followed by bin edges z_B(ntomo+1)"
  echo "       -t nbins number of theta bins, theta_min, theta_max"
  echo "       -i cross correlate bins i with j"
  echo "       -j cross correlate bins i with j"
  echo "       -c c-corr on? true/false"
  echo "       -l linear not log bins? true/false"
  echo "       -b which blind?"
  echo ""
  echo "DESCRIPTION: "
  echo "    The given mode corresponds to the 2pt stats or catalogue manipulation"
  echo "    step that you would like to run. Available modes are currently:"
  echo "      \"CREATETOMO\": Cut catalogues into tomographic bins and calculate and subtract c-term"
  echo ""            
  echo "      \"XI\": calculate xi+/- for tomo bin combination i j"
  echo ""  
  echo "      \"COMBINE\": combine the XI results from N and S for cross bin combination i j"
  echo ""           
  echo "      \"COSEBIS\": calculate En/Bn for tomo bin combination i j "
  echo ""           
  echo "      \"Pkk\": calculate cosmic shear Band powers for tomo bin combination i j "
  echo ""           
  echo "      \"GAMMAT\": calculate gamma_t and gamma_x for cross bin combination i j"
  echo ""           
  echo "      \"Pgk\": calculate GGL Band powers to cross bin combination i j"
  echo ""       
  echo "IMPORTANT DEPENDENCIES:"
  echo "    This script uses TreeCorr version 4.0.  Previous versions do not have linear binning"
  echo "    which is essential for COSEBIS"
  echo ""
  echo "EXAMPLES:"
  echo "    ./doall_calc2pt.sh -m \"CREATETOMO\""
  echo "        runs the CREATETOMO mode on the default data path and filters"
  echo ""
  echo ""
  echo "AUTHOR:"
  echo "    Catherine Heymans (heymans@roe.ac.uk)"
  echo ""
}

## Defaults
# Main Directory where the master catalogues are stored
MD=/disk09/KIDS/KIDSCOLLAB_V1.0.0/K1000_CATALOGUES_PATCH     
#Output Directory
OD=/disk09/KIDS/K1000_TWO_PT_STATS/   
# Catalogue Version numbery
LENSFIT_VER=v3             
# Analyse either North or South - and use the COMBINE mode 
# to combine the results.  Can be N, S, ALL     
PATCH=N
# Information about the tomographic bins
# Format:  ntomo, zb_edges (ntomo+ 1)
TOMOINFO="6 0.1 0.3 0.5 0.7 0.9 1.2 2.0"
# Information about the theta bins
# Format:  nbins, theta_min, theta_max
BININFO="9 0.5 300"
# Use a wrapper script to run over different bin 
# combinations - making it easier to run in parrallel
# This default correlates bin i=1 with bin j=2
# For GGL the first bin (i) is the lens bin, 
# and the second bin (j) is the source bin
IZBIN=1
JZBIN=2
# Do you want to apply a c-correction
CCORR=true
# Do you want to use linear binning
LINNOTLOG=false
# Which blind do you want to use?
BLIND=A


# Parse command line arguments
MODE=""

while getopts ":d:o:p:m:v:n:t:i:j:c:" opt; do
  case $opt in
    d)
      MD=$OPTARG
      ;;
    o) 
      OD=$OPTARG
      ;;
    p)
      PATCH=$OPTARG
      ;;
    m)
      MODE="$OPTARG"
      ;;
    v)
      LENSFIT_VER=$OPTARG
      ;;
    n)
      TOMOINFO=$OPTARG
      ;;
    t)
      BININFO=$OPTARG
      ;;
    i)
      IZBIN=$OPTARG
      ;;
    j)
      JZBIN=$OPTARG
      ;;
    c)
      CCORR=$OPTARG
      ;;
    l)
      LINNOTLOG=$OPTARG
      ;;
    b)
      BLIND=$OPTARG
      ;;
    
  esac
done

# If no MODE has been supplied then print usage and exit.
if [ -z "${MODE// }" ]; then
  echo "No mode was supplied. Please check the usage."
  printUsage
  exit 1
fi

##=================================================================
##
## Define some environment variables that are used in several modes
# If the file structures/names change, these will need editing

# STATDIR is the directory where all the results will go
STATDIR=${OD}/OUTSTATS
#Ensure all the directories that we want exist
mkdir -p $STATDIR/TOMOCATS
mkdir -p $STATDIR/XI

# And we're going to make some TMP files along the way that we'll want to easily delete so
USER=`whoami`
TMPDIR=/home/$USER/TMPDIR/2ptStats/
mkdir -p $TMPDIR

# The set-up below is for data
# We will need different file names for running through Flinc

# MASTERCAT is the main KiDS catalogue
MASTERCAT=${MD}/K1000_${PATCH}_9band_mask_BLINDED_${LENSFIT_VER}.cat

# TOMOCAT is the root name for the tomographic KiDS catalogue
# created by mode CREATETOMO

# To make life easier we will name our tomo bins with integer numbers 1,2,3,4,5,6
# Instead of the ZB boundaries
# We want to use these scripts for 2D aswell though, so we will preface
# with the total number of tomobins in the analysis

# Extract the number of tomobins from the TOMOINFO
TOMOINFOARR=($TOMOINFO)
NTOMO=${TOMOINFOARR[0]}

if [ $CCORR = "true" ]; then
  TOMOCAT=${OD}/TOMOCATS/K1000_${PATCH}_BLIND_${BLIND}_${LENSFIT_VER}_${NTOMO}Z
  C_RECORD=${OD}/TOMOCATS/c_terms_K1000_${PATCH}_BLIND_${BLIND}_${LENSFIT_VER}_${NTOMO}Z.asc
else
  TOMOCAT=${OD}/TOMOCATS/K1000_${PATCH}_BLIND_${BLIND}_${LENSFIT_VER}_NOCCORR_${NTOMO}Z
  C_RECORD=$TMPDIR/emptyfile  # just an empty file sent to the TMPDIR
fi

# Define the name for our output ascii file from Treecorr
BININFOARR=($BININFO)
outxi=$OD/XI/XI_K1000_${PATCH}_nbins_${BININFOARR[0]}_theta_${BININFOARR[1]}_${BININFOARR[2]}_zbins_${IZBIN}_${JZBIN}.asc

##=================================================================
##
## Now the description of modes for each catalogue processing task begins
## 
##=================================================================
##
#
## CREATETOMO\": Cut catalogues into tomographic bins"
#
for mode in ${MODE}
do
  if [ "$mode" = "CREATETOMO" ]; then
  
    echo "Starting mode CREATETOMO to cut catalogues into $NTOMO tomographic bins"

    if [ $PATCH = "ALL" ]; then
      { echo "MODE CREATETOMO only runs on PATCH N or S! Run MODE CREATETOMO -p N or -p S!"; exit 1; }
    fi

    # Check that the Master catalogue exists and exit if it doesn't 
    test -f ${MASTERCAT} || \
      { echo "Error: Master catalogue ${MASTERCAT} does not exist! Exiting!"; exit 1; }

    # Run through the NTOMO tomobins
    for ((i = 1 ; i <= $NTOMO ; i++)); do
        j=$(($i + 1))
        zmin=${TOMOINFOARR[$i]}
        zmax=${TOMOINFOARR[$j]}

        # script to select galaxies between zmin/zmax from ${MASTERCAT} and write out to ${TOMOCAT}_$i.cat
        # If CCORR is true, subtract off a c-term from the e1/e2 columns 
        #Also returns the correction applied to stdout which we send to $C_RECORD
        python create_tomocats.py $zmin $zmax ${MASTERCAT} ${TOMOCAT}_$i.fits $BLIND $CCORR

        # Check that the tomographic catalogues have been created and exit if they don't 
        test -f ${TOMOCAT}_$i.fits || \
        { echo "Error: Tomographic catalogue ${TOMOCAT}_$i.fits have not been created!"; exit 1; }
    done > $C_RECORD

    echo "Success: Leaving mode CREATETOMO"
  fi
done
##=================================================================
#
## \"XI\": calculate xi+/- for tomo bin combination i j"
#
for mode in ${MODE}
do
  if [ "$mode" = "XI" ]; then

    echo "Starting mode XI: calculate xi+/- for tomo bin combination $IZBIN $JZBIN with bins $BININFO"

    # Check that the tomographic catalogue exist and exit if they don't 
    test -f ${TOMOCAT}_$IZBIN.fits || \
      { echo "Error: Tomographic catalogue ${TOMOCAT}_$IZBIN.fits does not exist! Run MODE CREATETOMO!"; exit 1; }
    test -f ${TOMOCAT}_$JZBIN.fits || \
      { echo "Error: Tomographic catalogue ${TOMOCAT}_$JZBIN.fits does not exist! Run MODE CREATETOMO!"; exit 1; }

    # Run treecorr
    python calc_xi_w_treecorr.py $BININFO $LINNOTLOG ${TOMOCAT}_$IZBIN.fits ${TOMOCAT}_$JZBIN.fits $outxi

    # Did it work?
    test -f $outxi || \
      { echo "Error: Treecorr output $outxi was not created! !"; exit 1; }
    echo "Success: Leaving mode XI"

  fi
done

##=================================================================
#
#   \"COMBINE\": combine the results from N and S for cross bin combination i j"

for mode in ${MODE}
do
  if [ "$mode" = "COMBINE" ]; then

    echo "Starting mode COMBINE: to combine the N/S results for tomo bin \
          combination $IZBIN $JZBIN with bins $BININFO"

    # check do the files exist?
    tail=nbins_${BININFOARR[0]}_theta_${BININFOARR[1]}_${BININFOARR[2]}_zbins_${IZBIN}_${JZBIN}.asc
    outxiN=$OD/XI/XI_K1000_N_$tail
    outxiS=$OD/XI/XI_K1000_S_$tail

    test -f ${outxiN} || \
    { echo "Error: KiDS-N XI results $outxiN do not exist. Run MODE XI -p N!"; exit 1; } 
    test -f ${outxiS} || \
    { echo "Error: KiDS-S XI results $outxiS do not exist. Run MODE XI -p S!"; exit 1; } 

    # and now lets combine them using the fabulous awk
    # which Tilman will be most scathing about, but I love it nevertheless
    # first lets grab the header which we want to replicate

    head -1 < $outxiN > $TMPDIR/xi_header

    # paste the two catalogues together
    
    paste $outxiN $outxiS > $TMPDIR/xi_paste

    # time for awk where we use npairs to weight every other
    # column to get the average
    # $10 = npairs in KiDS-N,  $20 = npairs in KiDS-S
    # For the sigma_xi column I'm assuming ngals in N and S are similar and sum sigma_xi in quadrature
    # This isn't correct but we don't really use the sigma_xi column ($8 and $18) 
    # Finally give the sum of weights and the sum of npairs
    
    awk 'NR>1 {printf "%7.4e   %7.4e   %7.4e   %7.4e   %7.4e   %7.4e   %7.4e   %7.4e   %7.4e   %7.4e\n", 
                     ($1*$10 + $11*$20)/($10+$20), \
                     ($2*$10 + $12*$20)/($10+$20), \
                     ($3*$10 + $13*$20)/($10+$20), \
                     ($4*$10 + $14*$20)/($10+$20), \
                     ($5*$10 + $15*$20)/($10+$20), \
                     ($6*$10 + $16*$20)/($10+$20), \
                     ($7*$10 + $17*$20)/($10+$20), \
                     sqrt($8*$8 + $18*$18), $9+$19, $10+$20}' < $TMPDIR/xi_paste > $TMPDIR/xi_comb
    
    #finally put the header back

    outxi=$OD/XI/XI_K1000_ALL_$tail
    cat $TMPDIR/xi_header $TMPDIR/xi_comb > $outxi

    # Did it work?
    test -f $outxi || \
      { echo "Error: Combined Treecorr output $outxi was not created! !"; exit 1; }
    echo "Success: Leaving mode COMBINE"

  fi
done



##=================================================================
# To be written
#  echo ""           
#  echo "      \"COSEBIS\": calculate En/Bn for tomo bin combination i j "
#  echo ""           
#  echo "      \"Pkk\": calculate cosmic shear Band powers for tomo bin combination i j "
#  echo ""           
#  echo "      \"GAMMAT\": calculate gamma_t and gamma_x for cross bin combination i j"
#  echo ""           
#  echo "      \"Pgk\": calculate GGL Band powers to cross bin combination i j"