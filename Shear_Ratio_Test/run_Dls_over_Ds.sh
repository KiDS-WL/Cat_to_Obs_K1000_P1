#!/bin/bash

nlens=5                              # How many bins the lenses are divided into 
SOURCE_TYPE="K1000"                  # K1000 or MICE2_KV450
LENS_TYPE="GAMA_data"                # BOSS_data or MICE2_BOSS/MICE2_GAMA

                                     # ONLY USED FOR ANALYSIS WITH MICE
Realisation="Fid"                    # Which MICE realisation to use?
                                     # "Fid" is the single 343sqdeg realisation 
Mag_OnOff="on"                       # only used if SOURCE_TYPE is MICE2_KV450
Pz="Estimated"                       # only used if SOURCE_TYPE is MICE2_KV450

                                     # ONLY USED FOR ANALYSIS WITH K1000
SOMFLAGNAME="Fid"                    # The type of SOM used TO ESTIMATE n(Z)

nofz_shift=""             # Shift the nofz by the delta-z param?
                                      # Blank for no shift, "_nofzUp" for shift up, "_nofzDown" for shift down
                                      # "_nofzUp5sig" to specify 5sigma shift up (equivalent for Down).
                                      # "_nofzMix5sig" does incoherent (1,3,5 UP, 2,4 DOWN) 5 sigma shift.
                                      # nofz_shift SHOULD ONLY BE SET FOR K1000 SOURCES.

OL_Tag=""   #""   #"_OutlierPeaksInBins12" # If set, accesses nofzs saved in SOUCE_CATS/SOM_NofZ
                              # which have an Outlier (OL) excess injected at high redshift
                              # in the specified bin.
                          # SHOULD ONLY BE SET FOR K1000 SOURCES.
BLIND="C"                 # Which Blind to use


if [ "$SOURCE_TYPE" == "K1000" ]; then
    OUTDIR=Dls_over_Ds_data/SOURCE-${SOURCE_TYPE}_LENS-${LENS_TYPE}_Blind${BLIND}_SOM${SOMFLAGNAME}${OL_Tag}
elif [ "$SOURCE_TYPE" == "MICE2_KV450" ]; then

    if [ "$Realisation" == "Fid" ]; then
	OUTDIR=Dls_over_Ds_data/SOURCE-${SOURCE_TYPE}_LENS-${LENS_TYPE}_Pz${Pz}_mag${Mag_OnOff}
    elif [ "$Realisation" == "Octant" ]; then
	OUTDIR=Dls_over_Ds_data/SOURCE-${SOURCE_TYPE}_LENS-${LENS_TYPE}_Pz${Pz}_mag${Mag_OnOff}_Octant
    else
        echo "Realisation must be set to Fid or Octant. Currently it's set to $Realisation. EXITING."
        exit
    fi
	
fi
    
if [ ! -d ${OUTDIR} ]; then
    mkdir -p $OUTDIR
fi

    
JZBIN=1
for tomochar in 01t03 03t05 05t07 07t09 09t12
do
    if [ "$SOURCE_TYPE" == "K1000" ]; then
	# This is the old DIR-Estimated n(z)
	#source_nofz=/home/cech/public_html/KiDS/KiDS-VIKING-450/BLINDED_NZ/KiDS_2018-07-26_deepspecz_photoz_10th_BLIND_specweight_1000_4_ZB${tomochar}_blind${BLIND}_Nz.asc

	# This is the new SOM-Estimated n(z)
	if [ "$OL_Tag" == "" ]; then
	    source_nofz_DIR=/disk09/KIDS/K1000_TWO_PT_STATS/SOM_NofZ
	else
	    source_nofz_DIR=/home/bengib/KiDS1000_NullTests/Codes_4_KiDSTeam_Eyes/Shear_Ratio_Test/SOURCECATS/SOM_NofZ
	fi
	source_nofz=${source_nofz_DIR}/K1000_NS_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_svn_309c_2Dbins_v2_DIRcols_${SOMFLAGNAME}_blind${BLIND}_TOMO${JZBIN}_Nz${OL_Tag}.fits
	LENS_DIR=LENSCATS/${LENS_TYPE}/
	
    elif [ "$SOURCE_TYPE" == "MICE2_KV450" ]; then

	if [ "$Realisation" == "Fid" ]; then
	    source_nofz=/home/bengib/MICE2_Mocks/MICE2_KV450/nofz_CATS_Pz${Pz}/Nz_mag${Mag_OnOff}_120bins_${tomochar}.asc
	    LENS_DIR=LENSCATS/${LENS_TYPE}_mag${Mag_OnOff}/
	else
	    source_nofz=/home/bengib/MICE2_Mocks/MICE2_KV450/nofz_CATS_Pz${Pz}_Octant/Nz_mag${Mag_OnOff}_120bins_${tomochar}.asc
	    LENS_DIR=LENSCATS/${LENS_TYPE}_mag${Mag_OnOff}_Octant/
	fi
	
    fi
    
 
    for IZBIN in 1 2 3 4 5 
    do
	echo "----- Running Dls_over_Ds.py with Source bin $JZBIN, Lens bin $IZBIN -----"
	specz_fitsfile=${LENS_DIR}/lens_cat_${nlens}Z_${IZBIN}.fits
        outfile=${OUTDIR}/Dls_over_Ds_DIR_6Z_source_${JZBIN}_${nlens}Z_lens_${IZBIN}${nofz_shift}.asc
        python Dls_over_Ds.py $specz_fitsfile $source_nofz $nofz_shift > $outfile
    done
    JZBIN=$[$JZBIN +1]
done

