#!/bin/bash

JZBIN=0 
#for JZBIN in 1 2 3 4 5 
for tomochar in 01t03 03t05 05t07 07t09 09t12
do
    source_nofz=/home/cech/public_html/KiDS/KiDS-VIKING-450/BLINDED_NZ/KiDS_2018-07-26_deepspecz_photoz_10th_BLIND_specweight_1000_4_ZB${tomochar}_blindB_Nz.asc
    JZBIN=$[$JZBIN +1]
 
    for IZBIN in 1 2 3 4 5 
    do  

        specz_fitsfile=/disk09/KIDS/K1000_TWO_PT_STATS/GGLCATS/BOSS_data_5Z_${IZBIN}.fits
        outfile=Dls_over_Ds_data/Dls_over_Ds_DIR_6Z_source_${JZBIN}_5Z_lens_${IZBIN}.asc
        python Dls_over_Ds.py $specz_fitsfile $source_nofz > $outfile
    done    
done

