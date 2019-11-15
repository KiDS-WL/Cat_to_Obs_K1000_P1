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

# File inclusions:
# This file has the links to where your version of Python lives
# We expect treecorr to be installed as part of this
. ./progs.ini

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
  echo "       -d /path/to/source_catalogues"
  echo "       -g /path/to/lens_catalogues"
  echo "       -o /path/to/results"
  echo "       -p patch name N or S"
  echo "       -m list of modes"
  echo "       -v lensfit version"
  echo "       -n number of tomographic source bins, followed by bin edges z_B(ntomo+1)"
  echo "       -t number of theta bins, theta_min, theta_max"
  echo "       -a number of BP ell bins, ell_min, ell_max, apodisation"
  echo "       -e BP & COSEBIS theta_min, theta_max"
  echo "       -i cross correlate bins i with j - for GGL i is the lens bin"
  echo "       -j cross correlate bins i with j - for GGL j is the source bin"
  echo "       -c c-corr on? true/false"
  echo "       -l linear not log bins? true/false"
  echo "       -b which blind?"
  echo "       -u user catalogues"
  echo ""
  echo "DESCRIPTION: "
  echo "    The given mode corresponds to the 2pt stats or catalogue manipulation"
  echo "    step that you would like to run. Available modes are currently:"
  echo "      \"CREATETOMO\": cut catalogues into tomographic bins and calculate and subtract c-term"
  echo ""
  echo "      \"GT\": calculate gamma_t and gamma_x for cross bin combination i j"
  echo ""
  echo "      \"XI\": calculate xi+/- for tomo bin combination i j"
  echo ""
  echo "      \"Pgk\": calculate GGL bandpower to cross bin combination i j"
  echo ""
  echo "      \"Pkk\": calculate cosmic shear bandpower for tomo bin combination i j "
  echo ""
  echo "      \"COSEBIS\": calculate En/Bn for tomo bin combination i j "
  echo ""
  echo "      \"COMBINEGT\": combine the GT results from N and S for cross bin combination i j"
  echo ""
  echo "      \"COMBINEXI\": combine the XI results from N and S for cross bin combination i j"
  echo ""
  echo "IMPORTANT DEPENDENCIES:"
  echo "    This script uses TreeCorr version 4.0 to allow for linear or log binning"
  echo ""
  echo "EXAMPLES:"
  echo "    ./doall_calc2pt.sh -m \"CREATETOMO\""
  echo "        runs the CREATETOMO mode on the default data path and filters"
  echo ""
  echo "AUTHOR:"
  echo "    Catherine Heymans (heymans@roe.ac.uk)"
  echo ""
}

## Defaults

## Main Directory where the master KIDS catalogues are stored
MD=/disk09/KIDS/KIDSCOLLAB_V1.0.0/K1000_CATALOGUES_PATCH/

## Main Directory where the GGL catalogues are stored
MD2=/disk09/KIDS/K1000_TWO_PT_STATS/GGLCATS/

## Output Directory
## WARNING: Linc's temporary modification for test
OD=/disk09/KIDS/K1000_TWO_PT_STATS/
#OD=/disk05/calin/91_Data/Cat_to_Obs_K1000_P1/Calc_2pt_Stats/

## Catalogue Version number
LENSFIT_VER=v3

## Analyse either North or South - and use the COMBINE mode
## to combine the results.  Can be N, S, ALL
PATCH=N

## Information about the tomographic bins
## Format:  ntomo, zb_edges (ntomo+ 1)
#TOMOINFO_STR="6 0.1 0.3 0.5 0.7 0.9 1.2 2.0"
TOMOINFO_STR="5 0.1 0.3 0.5 0.7 0.9 1.2"

## Information about the GT & XI theta bins
## Format:  nbins, theta_min, theta_max
THETAINFO_STR="300 0.24428841736054135 403.49549216938652"
## This gives exact edges at 0.5 and 300 arcmin with 259 bins across that space.
## There are 29 bins below 0.5 & 12 bins beyond 300.

## Information about the BP ell bins
## Format:  nbins, ell_min, ell_max, do apodisation
ELLINFO_STR="8 100.0 1500.0 false"

## Information about the COSEBIS theta bins
## Format:  theta_min, theta_max
BP_COSEBIS_THETAINFO_STR="0.5 300"

## Use a wrapper script to run over different bin 
## combinations - making it easier to run in parrallel
## This default correlates bin i=1 with bin j=2
## For GGL the first bin (i) is the lens bin, 
## and the second bin (j) is the source bin
IZBIN=1
JZBIN=2

## Do you want to apply a c-correction
CCORR=true

## Do you want to use linear binning
LINNOTLOG=false

## Which blind do you want to use?
BLIND=A

## Do you want to define the input catalogue yourself with the -u 
## user defined catalogue option - if yes we need to set
USERCAT=false

## Parse command line arguments
MODE=""

while getopts ":d:g:o:p:m:v:n:t:a:e:i:j:c:l:b:u:" opt; do
  case $opt in
    d)
      MD=$OPTARG
      ;;
    g)
      MD2=$OPTARG
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
      TOMOINFO_STR=$OPTARG
      ;;
    t)
      THETAINFO_STR=$OPTARG
      ;;
    a)
      ELLINFO_STR=$OPTARG
      ;;
    e)
      BP_COSEBIS_THETAINFO_STR=$OPTARG
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
    u)
      USERCAT=$OPTARG
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
## If the file structures/names change, these will need editing

## Ensure all the directories that we want exist
mkdir -p ${OD}/TOMOCATS
mkdir -p ${OD}/OUTSTATS

## STATDIR is the directory where all the results will go
STATDIR=${OD}/OUTSTATS
mkdir -p ${STATDIR}/XI
mkdir -p ${STATDIR}/Pkk
mkdir -p ${STATDIR}/Pgk
mkdir -p ${STATDIR}/GT
mkdir -p ${STATDIR}/COSEBIS

## And we're going to make some TMP files along the way that we'll want to easily delete
## Depending on where you're running this script, this can be either in /home or /data 
## To be defined in progs.ini
mkdir -p $TMPDIR

## The default is to work with the main KiDS catalogue
## but you might want to work on mock data, in which case you can 
## fix the name of the mock catalogue with the -u option which sets MASTERCAT
if [ $USERCAT = "false" ]; then
  ## User catalogue has not been defined - use KIDS
  #Phase 1 catalogue
  catTag=K1000_${PATCH}_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_${LENSFIT_VER}
  MASTERCAT=${MD}/$catTag.cat
  #Phase 0 catalogue
  #MASTERCAT=${MD}/K1000_${PATCH}_9band_mask_BLINDED_${LENSFIT_VER}.cat
  #catTag=K1000_${PATCH}_BLIND_${BLIND}_${LENSFIT_VER}
else
  ## Set the filename of the output files to be the same as the name of the input fits catalogue
  MASTERCAT=${MD}/$USERCAT
  catTag=$USERCAT ## TODO
fi

## Tomographic bins
## Convert bin info strings into arrays
TOMOINFO=($TOMOINFO_STR)
NTOMO=${TOMOINFO[0]}
tomoPairTag="zbins_${IZBIN}_${JZBIN}"

## To make life easier we will name our tomo bins with integer numbers 1,2,3,4,5,6
## Instead of the ZB boundaries
## We want to use these scripts for 2D as well though, so we will preface
## with the total number of tomobins in the analysis
if [ $CCORR = "true" ]; then
  TOMOCAT=${OD}/TOMOCATS/${catTag}_${NTOMO}Z
  C_RECORD=${OD}/TOMOCATS/c_terms_${catTag}_${NTOMO}Z.asc
else
  TOMOCAT=${OD}/TOMOCATS/${catTag}_NOCCORR_${NTOMO}Z
  C_RECORD=$TMPDIR/emptyfile  # just an empty file sent to the TMPDIR
fi
## TOMOCAT is the root name for the catalogues created by mode CREATETOMO

## theta bins
THETAINFO=($THETAINFO_STR)
N_theta=${THETAINFO[0]}
angTag="nbins_${THETAINFO[0]}_theta_${THETAINFO[1]}_${THETAINFO[2]}"

## ell bins
ELLINFO=($ELLINFO_STR)
doApo=${ELLINFO[3]}
ellTag="nbins_${ELLINFO[0]}_Ell_${ELLINFO[1]}_${ELLINFO[2]}"

## theta range for BP & COSEBIs
BP_COSEBIS_THETAINFO=($BP_COSEBIS_THETAINFO_STR)

##=================================================================
##
## Now the description of modes for each catalogue processing task begins
## 
##=================================================================
##
##  \"CREATETOMO\": cut catalogues into tomographic bins"
##

for mode in ${MODE}
do
  if [ "${mode}" = "CREATETOMO" ]; then
  
    echo "Starting mode CREATETOMO to cut catalogues into $NTOMO tomographic bins"

    if [ ${PATCH} = "ALL" ]; then
      { echo "MODE CREATETOMO only runs on PATCH N or S! Run MODE CREATETOMO -p N or -p S!"; exit 1; }
    fi

    # Check that the Master catalogue exists and exit if it doesn't 
    test -f ${MASTERCAT} || \
      { echo "Error: Master catalogue ${MASTERCAT} does not exist! Exiting!"; exit 1; }

    # Run through the NTOMO tomobins
    for ((i = 1 ; i <= $NTOMO ; i++)); do
        j=$(($i + 1))
        zmin=${TOMOINFO[$i]}
        zmax=${TOMOINFO[$j]}

        # script to select galaxies between zmin/zmax from ${MASTERCAT} and write out to ${TOMOCAT}_$i.cat
        # If CCORR is true, subtract off a c-term from the e1/e2 columns 
        # Also returns the correction applied to stdout which we send to $C_RECORD
        $P_PYTHON create_tomocats.py $zmin $zmax ${MASTERCAT} ${TOMOCAT}_$i.fits $BLIND $CCORR

        # Check that the tomographic catalogues have been created and exit if they don't 
        test -f ${TOMOCAT}_$i.fits || \
        { echo "Error: Tomographic catalogue ${TOMOCAT}_$i.fits have not been created!"; exit 1; }
    done > $C_RECORD

    echo "Success: Leaving mode CREATETOMO"
  fi
done

##==============================================================================
##
##  \"GT\": calculate gamma_t and gamma_x for cross bin combination i j"
##  We want to use this mode for both the shear-ratio test and the main 3x2pt
##  Analysis so it needs to be flexible
##

for mode in ${MODE}
do
  if [ "${mode}" = "GT" ]; then

    echo "Starting mode ${mode} to calculate gamma_t for lens bin ${IZBIN} + source bin ${JZBIN} \
          with a total number of tomo bins $THETAINFO_STR"

    if [ "${PATCH}" = "N" ]; then
       GGL_ID="BOSS"
    elif [ "${PATCH}" = "S" ]; then
       GGL_ID="2dFLenS"
    else
      { echo "MODE ${mode} only runs on PATCH N or S! Run MODE CREATETOMO -p N or -p S!"; exit 1; }
    fi
    
    lensCat="$MD2/${GGL_ID}_data_z$IZBIN.fits"
    randCat="$MD2/${GGL_ID}_random_z$IZBIN.fits"
    srcCat="${TOMOCAT}_$JZBIN.fits"

    # check does the correct lens/source/random files exist?
    test -f ${lensCat} || \
      { echo "Error: Lens catalogue ${lensCat} does not exist."; exit 1; } 
    test -f ${randCat} || \
      { echo "Error: Random catalogue ${randCat} does not exist."; exit 1; } 
    test -f ${srcCat} || \
      { echo "Error: Tomographic catalogue ${srcCat} does not exist! Run MODE CREATETOMO!"; exit 1; }

    # where do we want to write the output to?
    outPath=${STATDIR}/GT/GT_${catTag}_${angTag}_${tomoPairTag}.asc

    # Run treecorr - using the Mandelbaum estimator that subtracts off the random signal
    $P_PYTHON calc_gt_w_treecorr.py $THETAINFO_STR $LINNOTLOG ${lensCat} ${randCat} ${srcCat} ${outPath}

    # Did it work?
    test -f ${outPath} || \
      { echo "Error: Treecorr output ${outPath} was not created! !"; exit 1; }
    echo "Success: Leaving mode ${mode}"
  fi
done

##==============================================================================
##
##  \"XI\": calculate xi+/- for tomo bin combination i j"
##

for mode in ${MODE}
do
  if [ "${mode}" = "XI" ]; then

    echo "Starting mode ${mode}: calculate xi+/- for tomo bin combination $IZBIN $JZBIN with bins $THETAINFO_STR"

    srcCat1="${TOMOCAT}_$IZBIN.fits"
    srcCat2="${TOMOCAT}_$JZBIN.fits"

    # Check that the tomographic catalogue exist and exit if they don't 
    test -f ${srcCat1} || \
      { echo "Error: Tomographic catalogue ${srcCat1} does not exist! For KiDS run MODE CREATETOMO!"; exit 1; }
    test -f ${srcCat2} || \
      { echo "Error: Tomographic catalogue ${srcCat2} does not exist! For KiDS run MODE CREATETOMO!"; exit 1; }

    # Define the name for our output ascii files from Treecorr etc
    outPath=${STATDIR}/XI/XI_${catTag}_${angTag}_${tomoPairTag}.asc

    # Run treecorr
    $P_PYTHON calc_xi_w_treecorr.py $THETAINFO_STR $LINNOTLOG ${srcCat1} ${srcCat2} ${outPath}

    # Did it work?
    test -f ${outPath} || \
      { echo "Error: Treecorr output ${outPath} was not created! !"; exit 1; }
    echo "Success: Leaving mode ${mode}"
  fi
done

##==============================================================================
##
##  \"Pgk\": calculate GGL bandpower to cross bin combination i j"
##

for mode in ${MODE}
do
  if [ "${mode}" = "Pgk" ]; then

    echo "Starting mode ${mode}: to calculate GGL bandpower for lens bin ${IZBIN} + source bin ${JZBIN} \
          with bins ${ELLINFO_STR} & theta range ${BP_COSEBIS_THETAINFO_STR}"

    # Check whether the correct xi files exist?
    InputFileIdentifier=${catTag}_${angTag}_${tomoPairTag}
    treePath=${STATDIR}/XI/XI_${InputFileIdentifier}.asc

    test -f ${treePath} || \
    { echo "Error: KiDS-${PATCH} GT results ${treePath} do not exist. Either Run MODE GT (N/S) or COMBINE (ALL)!"; exit 1; } 

    ## Define correlation type (1: ee; 2: ne; 3: nn)
    corrType=2

    # Define the workplace where crazy file handling happens
    FolderName=${STATDIR}/Pgk

    # The files need to have this naming convention:  xi2bandpow_input_${InputFileIdentifier}.dat
    # They also need to have only 3 columns with no other lines or comments:
    # theta [arcmin]    [gamma_t or xi_+]    [gamma_x or xi_-]

    # Let's use awk to convert the Treecorr output into the expected format and remove the header.
    # We will use the meanR as the radius to pass through to the bandpower code
    # This choise is not too important though the bins are finely binned
    # Treecorr: #   R_nom       meanR       meanlogR       xip          xim         xip_im      xim_im      sigma_xi      weight       npairs

    if [ "${doApo}" = "true" ]; then
        ## apoWidth = log width of apodisation window
        ## If nominal theta range (BP_COSEBIS_THETAINFO_STR) is theta_min to theta_max,
        ## the real theta range that input file should provide is real_min to real_max, 
        ## where real_min = theta_min / exp(apoWidth/2),
        ##       real_max = theta_max * exp(apoWidth/2).
        awk '(NR>1){print $2, $4-$5, $10-$11}' < ${treePath} > ${FolderName}/xi2bandpow_input_${InputFileIdentifier}.dat
        apoWidth=0.5
        N_theta_BP=300 # WARNING Temporary
    else
        awk '(NR>30 && NR<=289){print $2, $4-$5, $10-$11}' < ${treePath} > ${FolderName}/xi2bandpow_input_${InputFileIdentifier}.dat
        apoWidth=-1
        N_theta_BP=259 # WARNING Temporary
    fi

    # The output file is called xi2bandpow_output_${OutputFileIdentifier}.dat 
    # It has 3 columns for GGL: ell, bandpow, err
    # And 5 columns for cosmic shear: ell, l*2(P_E/2pi), err, l*2(P_B/2pi), err
    # ell is the log-central value
    OutputFileIdentifier=${catTag}_${ellTag}_${tomoPairTag}

    # These are the options for inputs for the c program xi2bandpow.c:
    # 1: <working directory>
    # 2: <input file identifier>
    # 3: <output file identifier>
    # 4: <number of input angular bins>
    # 5: <min input separation to use in conversion [arcmin] (xi_+ or gamma_+)>
    # 6: <max input separation to use in conversion [arcmin] (xi_+ or gamma_+)>
    # 7: <min input separation to use in conversion [arcmin] (xi_- or gamma_x)>
    # 8: <max input separation to use in conversion [arcmin] (xi_- or gamma_x)>
    # 9: <number of output ell bins>
    # 10: <min output ell>
    # 11: <max output ell>
    # 12: <correlation type (1: ee; 2: ne; 3: nn)>
    # 13: <log width of apodisation window [total width of apodised range is tmax/tmin=exp(width) in arcmin; <0 for no apodisation]>
    # now run the program (location is stored in progs.ini)
    $P_XI2BANDPOW ${FolderName} ${InputFileIdentifier} ${OutputFileIdentifier} \
                  ${N_theta_BP} ${BP_COSEBIS_THETAINFO[0]} ${BP_COSEBIS_THETAINFO[1]} ${BP_COSEBIS_THETAINFO[0]} ${BP_COSEBIS_THETAINFO[1]} \
                  ${ELLINFO[0]} ${ELLINFO[1]} ${ELLINFO[2]} ${corrType} ${apoWidth}

    rm -f "${FolderName}/xi2bandpow_input_*"
    rm -f "${FolderName}/xi2bandpow_kernels_*"
    outPath=${FolderName}/xi2bandpow_output_${OutputFileIdentifier}.dat

    # Did it work?
    test -f ${outPath} || \
      { echo "Error: bandpower output ${outPath} was not created! !"; exit 1; }
    echo "Success: Leaving mode ${mode}"

  fi
done

##==============================================================================
##
##  \"Pkk\": calculate cosmic shear bandpower for tomo bin combination i j "

for mode in ${MODE}
do
  if [ "${mode}" = "Pkk" ]; then

    echo "Starting mode ${mode}: to calculate cosmic shear Band powers for tomo bin combination \
          ${IZBIN} ${JZBIN} with bins ${ELLINFO_STR} & theta range ${BP_COSEBIS_THETAINFO_STR}"

    ## Check whether the correct xi files exist?
    InputFileIdentifier=${catTag}_${angTag}_${tomoPairTag}
    treePath=${STATDIR}/XI/XI_${InputFileIdentifier}.asc

    test -f ${treePath} || \
    { echo "Error: KiDS-${PATCH} XI results ${treePath} do not exist. Either Run MODE XI (N/S) or COMBINE (ALL)!"; exit 1; } 

    ## Define correlation type (1: ee; 2: ne; 3: nn)
    corrType=1

    # Define the workplace where crazy file handling happens
    FolderName=${STATDIR}/Pkk

    # The files need to have this naming convention:  xi2bandpow_input_${InputFileIdentifier}.dat
    # They also need to have only 3 columns with no other lines or comments:
    # theta [arcmin]    [gamma_t or xi_+]    [gamma_x or xi_-]

    # Let's use awk to convert the Treecorr output into the expected format and remove the header.
    # We will use the meanR as the radius to pass through to the bandpower code
    # This choise is not too important though the bins are finely binned
    # Treecorr: #   R_nom       meanR       meanlogR       xip          xim         xip_im      xim_im      sigma_xi      weight       npairs

    if [ "${doApo}" = "true" ]; then
        ## apoWidth = log width of apodisation window
        ## If nominal theta range (BP_COSEBIS_THETAINFO_STR) is theta_min to theta_max,
        ## the real theta range that input file should provide is real_min to real_max, 
        ## where real_min = theta_min / exp(apoWidth/2),
        ##       real_max = theta_max * exp(apoWidth/2).
        awk '(NR>2){print $2, $4, $5}' < ${treePath} > ${FolderName}/xi2bandpow_input_${InputFileIdentifier}.dat
        apoWidth=0.5
        N_theta_BP=300 # WARNING Temporary
    else
        awk '(NR>31 && NR<=290){print $2, $4, $5}' < ${treePath} > ${FolderName}/xi2bandpow_input_${InputFileIdentifier}.dat
        apoWidth=-1
        N_theta_BP=259 # WARNING Temporary
    fi

    # The output file is called xi2bandpow_output_${OutputFileIdentifier}.dat 
    # It has 3 columns for GGL: ell, bandpow, err
    # And 5 columns for cosmic shear: ell, l*2(P_E/2pi), err, l*2(P_B/2pi), err
    # ell is the log-central value
    OutputFileIdentifier=${catTag}_${ellTag}_${tomoPairTag}

    # These are the options for inputs for the c program xi2bandpow.c:
    # 1: <working directory>
    # 2: <input file identifier>
    # 3: <output file identifier>
    # 4: <number of input angular bins>
    # 5: <min input separation to use in conversion [arcmin] (xi_+ or gamma_+)>
    # 6: <max input separation to use in conversion [arcmin] (xi_+ or gamma_+)>
    # 7: <min input separation to use in conversion [arcmin] (xi_- or gamma_x)>
    # 8: <max input separation to use in conversion [arcmin] (xi_- or gamma_x)>
    # 9: <number of output ell bins>
    # 10: <min output ell>
    # 11: <max output ell>
    # 12: <correlation type (1: ee; 2: ne; 3: nn)>
    # 13: <log width of apodisation window [total width of apodised range is tmax/tmin=exp(width) in arcmin; <0 for no apodisation]>
    # now run the program (location is stored in progs.ini)
    $P_XI2BANDPOW ${FolderName} ${InputFileIdentifier} ${OutputFileIdentifier} \
                  ${N_theta_BP} ${BP_COSEBIS_THETAINFO[0]} ${BP_COSEBIS_THETAINFO[1]} ${BP_COSEBIS_THETAINFO[0]} ${BP_COSEBIS_THETAINFO[1]} \
                  ${ELLINFO[0]} ${ELLINFO[1]} ${ELLINFO[2]} ${corrType} ${apoWidth}

    rm -f "${FolderName}/xi2bandpow_input_*"
    rm -f "${FolderName}/xi2bandpow_kernels_*"
    outPath=${FolderName}/xi2bandpow_output_${OutputFileIdentifier}.dat

    # Did it work?
    test -f ${outPath} || \
      { echo "Error: bandpower output ${outPath} was not created! !"; exit 1; }
    echo "Success: Leaving mode ${mode}"

  fi
done

##==============================================================================
##
##  \"COSEBIS\": calculate COSEBIs"
##  The angular ranges that we can probe are limited by the pre-computed tables
##  in src/cosebis/TLogsRootsAndNorms/
##

for mode in ${MODE}
do
    if [ "${mode}" = "COSEBIS" ]; then

    echo "Starting mode ${mode}: to calculate COSEBIS for bin combination \
    ${IZBIN} ${JZBIN} on theta range ${BP_COSEBIS_THETAINFO_STR}"

    # check does the correct input xi file exist?
    treePath=${STATDIR}/XI/XI_${catTag}_${angTag}_${tomoPairTag}.asc

    test -f ${treePath} || \
    { echo "Error: KiDS-${PATCH} XI results ${treePath} do not exist. Either Run MODE XI (N/S) or COMBINE (ALL)!"; exit 1; }

    # check that the pre-computed COSEBIS tables exist
    SRCLOC=../src/cosebis
    normfile=$SRCLOC/TLogsRootsAndNorms/Normalization_${BP_COSEBIS_THETAINFO[0]}-${BP_COSEBIS_THETAINFO[1]}.table
    rootfile=$SRCLOC/TLogsRootsAndNorms/Root_${BP_COSEBIS_THETAINFO[0]}-${BP_COSEBIS_THETAINFO[1]}.table

    test -f ${normfile} || \
    { echo "Error: COSEBIS pre-computed table $normfile is missing. Download from gitrepo!"; exit 1; }

    test -f ${rootfile} || \
    { echo "Error: COSEBIS pre-computed table $rootfile is missing. Download from gitrepo!"; exit 1; }

    # have we run linear or log binning for the 2pt correlation function?
    if [ "$LINNOTLOG" = "false" ]; then 
      binning='log'
    else
      binning='lin'
    fi

    # where do we want to write the output to?
    filetail=COSEBIS_${catTag}_theta_${BP_COSEBIS_THETAINFO[0]}_${BP_COSEBIS_THETAINFO[1]}_${tomoPairTag}
    outcosebis=${STATDIR}/COSEBIS/

    # Now Integrate output from treecorr with COSEBIS filter functions
    # -i = input file
    # -t = treecorr output theta_col - the first column is zero so -t 1 uses the meanR from Treecorr
    # -p = treecorr output xip_col
    # -m = treecorr output xim_col
    # --cfoldername = output directory
    # -o = filename (outputs En_filename.ascii and Bn_filename.ascii)
    # -b = binning "log" or "lin"
    # -n = number of COSEBIS modes
    # -s = COSEBIS minimum theta
    # -l = COSEBIS maximum theta
    # location of the required pre-compution tables
    # --tfoldername = Tplus_minus    # computes/saves the first time it is run for a fixed theta min/max
    # --norm = TLogsRootsAndNorms/Normalization_${tmin}-${tmax}.table
    # --root = TLogsRootsAndNorms/Root_${tmin}-${tmax}.table

    $P_PYTHON ../src/cosebis/run_measure_cosebis_cats2stats.py -i $xifile -t 1 -p 3 -m 4 \
            --cfoldername $outcosebis -o $filetail -b $binning -n 5 -s ${BP_COSEBIS_THETAINFO[0]} \
            -l ${BP_COSEBIS_THETAINFO[1]} --tfoldername $SRCLOC/Tplus_minus \
            --norm $normfile --root $rootfile

    # I am expecting this to have produced two files called
    Enfile=$outcosebis/En_${filetail}.asc
    Bnfile=$outcosebis/Bn_${filetail}.asc

    # Did it work?
    test -f $Enfile || \
    { echo "Error: COSEBIS measurement $Enfile was not created! !"; exit 1; }
    test -f $Bnfile || \
    { echo "Error: COSEBIS measurement $Bnfile was not created! !"; exit 1; }

    echo "Success: Leaving mode ${mode}"

    fi
done

##==============================================================================
##
##  \"COMBINEGT\": combine the GT results from N and S for cross bin combination i j"
##

for mode in ${MODE}
do
  if [ "${mode}" = "COMBINEGT" ]; then

    echo "Starting mode COMBINE: to combine the N/S results for tomo bin \
          combination $IZBIN $JZBIN with bins $THETAINFO_STR"

    # check do the files exist?
    inFileN=${STATDIR}/GT/GT_K1000_N_${angTag}_${tomoPairTag}.asc
    inFileS=${STATDIR}/GT/GT_K1000_S_${angTag}_${tomoPairTag}.asc

    test -f ${inFileN} || \
    { echo "Error: KiDS-N GT results ${inFileN} do not exist. Run MODE GT -p N!"; exit 1; } 
    test -f ${inFileS} || \
    { echo "Error: KiDS-S GT results ${inFileS} do not exist. Run MODE GT -p S!"; exit 1; } 

    # and now lets combine them using the fabulous awk
    # which Tilman will be most scathing about, but I love it nevertheless
    # first lets grab the header which we want to replicate

    head -1 < ${inFileN} > $TMPDIR/xi_header

    # paste the two catalogues together
    
    paste ${inFileN} ${inFileS} > $TMPDIR/xi_paste

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

    outPath=${STATDIR}/XI/XI_K1000_ALL_${angTag}_${tomoPairTag}.asc
    cat $TMPDIR/xi_header $TMPDIR/xi_comb > ${outPath}

    # Did it work?
    test -f ${outPath} || \
      { echo "Error: Combined Treecorr output ${outPath} was not created! !"; exit 1; }
    echo "Success: Leaving mode COMBINE"

  fi
done

##==============================================================================
##
##  \"COMBINEXI\": combine the XI results from N and S for cross bin combination i j"
##

for mode in ${MODE}
do
  if [ "${mode}" = "COMBINEXI" ]; then

    echo "Starting mode COMBINE: to combine the N/S results for tomo bin \
          combination $IZBIN $JZBIN with bins $THETAINFO_STR"

    # check do the files exist?
    #Phase 0
    #inFileN=${STATDIR}/XI/XI_K1000_N_BLIND_${BLIND}_${LENSFIT_VER}_${angTag}_${tomoPairTag}.asc
    #inFileS=${STATDIR}/XI/XI_K1000_S_BLIND_${BLIND}_${LENSFIT_VER}_${angTag}_${tomoPairTag}.asc

    #Phase 1
    inFileN=$STATDIR/XI/XI_K1000_N_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_${LENSFIT_VER}_${angTag}_${tomoPairTag}.asc
    inFileS=$STATDIR/XI/XI_K1000_S_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_${LENSFIT_VER}_${angTag}_${tomoPairTag}.asc

    test -f ${inFileN} || \
    { echo "Error: KiDS-N XI results ${inFileN} do not exist. Run MODE XI -p N!"; exit 1; } 
    test -f ${inFileS} || \
    { echo "Error: KiDS-S XI results ${inFileS} do not exist. Run MODE XI -p S!"; exit 1; } 

    # and now lets combine them using the fabulous awk
    # which Tilman will be most scathing about, but I love it nevertheless
    # first lets grab the header which we want to replicate

    # old versions of treecorr
    #head -1 < ${inFileN} > $TMPDIR/xi_header
    # latest version of treecorr (3.7)
    head -2 < ${inFileN} > $TMPDIR/xi_header

    # paste the two catalogues together
    paste ${inFileN} ${inFileS} > $TMPDIR/xi_paste

    # time for awk where we use npairs to weight every other
    # column to get the average
    # $11 = npairs in KiDS-N,  $22 = npairs in KiDS-S
    # For the sigma_xi column I'm assuming ngals in N and S are similar and sum sigma_xi in quadrature
    # This isn't correct but we don't really use the sigma_xi column ($8 and $18) 
    # Finally give the sum of weights and the sum of npairs
    
    awk 'NR>2 {printf "%7.4e   %7.4e   %7.4e   %7.4e   %7.4e   %7.4e   %7.4e   %7.4e   %7.4e   %7.4e %7.4e\n", 
                     ($1*$11 + $12*$22)/($11+$22), \
                     ($2*$11 + $13*$22)/($11+$22), \
                     ($3*$11 + $14*$22)/($11+$22), \
                     ($4*$11 + $15*$22)/($11+$22), \
                     ($5*$11 + $16*$22)/($11+$22), \
                     ($6*$11 + $17*$22)/($11+$22), \
                     ($7*$11 + $18*$22)/($11+$22), \
                     sqrt($8*$8 + $19*$19), sqrt($9*$9 + $20*$20), $10+$21, $11+$22}' < $TMPDIR/xi_paste > $TMPDIR/xi_comb

    #finally put the header back
    #Phase 0
    #outPath=${STATDIR}/XI/XI_K1000_ALL_${angTag}_${tomoPairTag}.asc
    #Phase 1
    outPath=$STATDIR/XI/XI_K1000_ALL_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_${LENSFIT_VER}_${angTag}_${tomoPairTag}.asc
    cat $TMPDIR/xi_header $TMPDIR/xi_comb > ${outPath}

    # Did it work?
    test -f ${outPath} || \
      { echo "Error: Combined Treecorr output ${outPath} was not created! !"; exit 1; }
    echo "Success: Leaving mode COMBINE"

  fi
done
