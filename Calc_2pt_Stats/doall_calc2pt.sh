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
  echo "       -s patch name N or S"
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
  echo "      \"COSEBIS\": calculate En/Bn for tomo bin combination i j "
  echo ""           
  echo "      \"Pkk\": calculate cosmic shear Band powers for tomo bin combination i j "
  echo ""           
  echo "      \"GAMMAT\": calculate gamma_t and gamma_x for cross bin combination i j"
  echo ""           
  echo "      \"Pgk\": calculate GGL Band powers to cross bin combination i j"
  echo ""           
  echo "      \"COMBINE\": combine the results from N and S for cross bin combination i j"
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
# to combine the results     
PATCH=N
# Information about the tomographic bins
# Format:  ntomo, zb_edges (ntomo+ 1)
TOMOINFO="6 0.1 0.3 0.5 0.7 0.9 1.2 2.0"
# Information about the theta bins
# Format:  nbins, theta_min, theta_max
BININFO="7 0.5 300"
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

while getopts ":d:o:s:m:v:n:t:i:j:c:" opt; do
  case $opt in
    d)
      MD=$OPTARG
      ;;
    o) 
      OD=$OPTARG
      ;;
    s)
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

# MASTERCAT is the main KiDS catalogue
#MASTERCAT=${MD}/K1000_${PATCH}_9band_mask_BLINDED_${LENSFIT_VER}.cat
MASTERCAT=${MD}/blinded_KIDS_9p0_m28p2.tmp

# TOMOCAT is the root name for the tomographic KiDS catalogue
# created by mode CREATETOMO

# To make life easier we will name our tomo bins with integer numbers 1,2,3,4,5
# Instead of the ZB boundaries
# We want to use these scripts for 2D aswell though, so we will preface
# with the total number of tomobins in the analysis

# Extract the number of tomobins from the TOMOINFO
TOMOINFOARR=($TOMOINFO)
NTOMO=${TOMOINFOARR[0]}

if [ $CCORR = "true" ]; then
  TOMOCAT=${OD}/TOMOCATS/K1000_${PATCH}_BLIND_${BLIND}_${LENSFIT_VER}_${NTOMO}Z
else
  TOMOCAT=${OD}/TOMOCATS/K1000_${PATCH}_BLIND_${BLIND}_${LENSFIT_VER}_NOCCORR_${NTOMO}Z
fi

# STATDIR is the directory where all the results will go
STATDIR=${OD}/OUTSTATS

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
        python create_tomocats.py $zmin $zmax ${MASTERCAT} ${TOMOCAT}_$i.cat $BLIND CCORR

        # Check that the tomographic catalogues have been created and exit if they don't 
        test -f ${TOMOCAT}_$i.cat || \
        { echo "Error: Tomographic catalogue ${TOMOCAT}_$i.cat have not been created!"; exit 1; }
    done

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
    test -f ${TOMOCAT}_$IZBIN.cat || \
      { echo "Error: Tomographic catalogue ${TOMOCAT}_$IZBIN.cat does not exist! Run MODE CREATETOMO!"; exit 1; }
    test -f ${TOMOCAT}_$JZBIN.cat || \
      { echo "Error: Tomographic catalogue ${TOMOCAT}_$JZBIN.cat does not exist! Run MODE CREATETOMO!"; exit 1; }

    # Set the name for our output ascii file from Treecorr
    BININFOARR=($BININFO)
    outxi=$OD/XI/XI_nbins_${BININFOARR[0]}_theta_${BININFOARR[1]}_${BININFOARR[2]}_zbins_${IZBIN}_${JZBIN}.asc

    # Run treecorr
    python calc_xi_w_treecorr.py $BININFO $LINNOTLOG ${TOMOCAT}_$IZBIN.cat ${TOMOCAT}_$JZBIN.cat $outxi

    # Did it work?
    test -f $outxi || \
      { echo "Error: Treecorr output $outxi was not created! !"; exit 1; }
    echo "Success: Leaving mode XI"

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
#  echo ""           
#  echo "      \"COMBINE\": combine the results from N and S for cross bin combination i j"