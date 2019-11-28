#!/bin/bash

# 20/11/2019, B. Giblin, Edinburgh
# script to run create_tomocats.py with bins specified below.
# NOTE: haven't yet changed Shear_ratio_wspin_test.py to read the source
# and lens catalogues saved by this script. They still read those previously
# saved in /disk09/KIDS/...


LENSCAT=/disk09/KIDS/K1000_TWO_PT_STATS/GGLCATS/BOSS_data_z1.fits
LENS_LO_EDGES=(0.2 0.3 0.4)
LENS_HI_EDGES=(0.3 0.4 0.5)  # Made a separate array for upper edges of bins
# since bash is chronically bad adding floats

for i in `seq 0 2`; do
    l="${LENS_LO_EDGES[i]}"
    h="${LENS_HI_EDGES[i]}"
    #python create_tomocats_GGL.py $l $h $LENSCAT Output/lens_tomocat_${l}_${h}.fits
done

SOURCECAT=/disk09/KIDS/KIDSCOLLAB_V1.0.0/K1000_CATALOGUES_PATCH/K1000_N_Phase0_9band_mask_svn_309b.cat
SOURCE_LO_EDGES=(0.1 0.3 0.5 0.7 0.9)
SOURCE_HI_EDGES=(0.3 0.5 0.7 0.9 1.2)
for i in `seq 0 5`; do
    l="${SOURCE_LO_EDGES[i]}"
    h="${SOURCE_HI_EDGES[i]}"
    python create_tomocats_GGL.py $l $h $SOURCECAT Output/source_tomocat_${l}_${h}.fits
done
