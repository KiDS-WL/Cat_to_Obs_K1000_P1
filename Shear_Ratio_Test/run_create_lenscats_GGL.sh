#!/bin/bash

# 20/11/2019, B. Giblin, Edinburgh
# script to run create_tomocats.py with bins specified below.



# INPUT PARAMATERS
# -----------------------------------------------------------------------
# Set the input catalogue based on the data type (BOSS-data/BOSS-randoms OR K1000/MICE mocks)

LENS_TYPE="GAMA_random"            # "BOSS_data", "BOSS_random",
                                  # "MICE2_BOSS", "MICE2_BOSS_random",
Realisation="Fid"                 # Which MICE realisation to use?
                                  # "Fid" is the single 343sqdeg realisation
                                  # "Octant" is the octant-sized realisation 
Mag_OnOff="on"                   # Only used for MICE.
# -----------------------------------------------------------------------

if [ "$LENS_TYPE" == "BOSS_data" ] || [ "$LENS_TYPE" == "BOSS_random" ]; then
    LENSCAT=/disk09/KIDS/K1000_TWO_PT_STATS/GGLCATS/${LENS_TYPE}_z1.fits
    LENS_OUTDIR=LENSCATS/${LENS_TYPE}/

    # Change the lens binning depending on the LENS_TYPE
    LENS_LO_EDGES=(0.2 0.3 0.4 0.5 0.6)
    LENS_HI_EDGES=(0.3 0.4 0.5 0.6 0.7)  # Made a separate array for upper edges of bins
                                         # since bash is chronically bad adding floats 

elif [ "$LENS_TYPE" == "GAMA_data" ] || [ "$LENS_TYPE" == "GAMA_random" ]; then
    LENS_LO_EDGES=(0.0 0.1 0.2 0.3 0.4)
    LENS_HI_EDGES=(0.1 0.2 0.3 0.4 0.5)
    LENS_OUTDIR=LENSCATS/${LENS_TYPE}/
    if [ "$LENS_TYPE" == "GAMA_data" ]; then
	LENSCAT=LENSCATS/${LENS_TYPE}/kids1000N_lenses_gama.dat
    elif [ "$LENS_TYPE" == "GAMA_random" ]; then
	LENSCAT=LENSCATS/${LENS_TYPE}/kids1000N_randlenses_gama_short.dat
    fi
    
    
elif [[ "$LENS_TYPE" == *"MICE2"* ]]; then
    LENS_TYPE_prefix=${LENS_TYPE%*"_random"*}   # Remove random tag if there.
    # Randomisation applied in create_lenscats_GGL.py

    # Change the lens binning depending on the LENS_TYPE
    if [[ "$LENS_TYPE" == *"BOSS"* ]]; then
	LENS_LO_EDGES=(0.2 0.3 0.4 0.5 0.6)
	LENS_HI_EDGES=(0.3 0.4 0.5 0.6 0.7)
    elif [[ "$LENS_TYPE" == *"GAMA"* ]]; then
	LENS_LO_EDGES=(0.0 0.1 0.2 0.3 0.4)
        LENS_HI_EDGES=(0.1 0.2 0.3 0.4 0.5)
    fi


    # Which MICE realisation to make the Lens Cat for? Fiducial (343sqdeg) or Octant-sized?
    if [ "$Realisation" == "Fid" ]; then
	LENSCAT=/home/bengib/MICE2_Mocks/MICE2_KV450/${LENS_TYPE_prefix}_magnification_${Mag_OnOff}_small.fits
	LENS_OUTDIR=LENSCATS/${LENS_TYPE}_mag${Mag_OnOff}/
    elif [ "$Realisation" == "Octant" ]; then
	LENSCAT=/disk09/KIDS/MICE_BOSS_KV_MOCKS/${LENS_TYPE_prefix}_magnification_${Mag_OnOff}.fits
        LENS_OUTDIR=LENSCATS/${LENS_TYPE}_mag${Mag_OnOff}_Octant/
    else
	echo "Realisation must be set to Fid or Octant. Currently it's set to $Realisation. EXITING."
	exit
    fi

    
fi


if [ ! -d "$LENS_OUTDIR" ]; then
    mkdir -p $LENS_OUTDIR
fi

# Used to use this binning for all Lens types. No longer.
#LENS_LO_EDGES=(0.2 0.3 0.4)
#LENS_HI_EDGES=(0.3 0.4 0.5)  # Made a separate array for upper edges of bins
                             # since bash is chronically bad adding floats

nzlens=${#LENS_LO_EDGES[@]} # number lens bins
for i in `seq 1 $nzlens`; do
    idx=$((i-1))
    l="${LENS_LO_EDGES[idx]}"
    h="${LENS_HI_EDGES[idx]}"
    python create_lenscats_GGL.py $l $h $LENSCAT $LENS_OUTDIR/lens_cat_${nzlens}Z_${i}.fits
done



