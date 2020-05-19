#!/bin/bash
#name_in = XI_K1000_ALL_BLIND_A_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_svn_309c_2Dbins_v2_goldclasses_Flag_SOM_Fid_nbins_326_theta_0.37895134266193781_395.82918204307509_zbins_
name_in=XI_K1000_ALL_BLIND_A_V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_svn_309c_2Dbins_v2_goldclasses_Flag_SOM_Fid_nbins_4000_theta_0.5_300.0_zbins_
name_out=XI_V1_4000
for r1 in `seq 1 5`; do
	for r2 in `seq ${r1} 5`; do
		cp ${name_in}${r1}_${r2}.asc ${name_out}_nBins_5_Bin${r1}_Bin${r2}.ascii
	done
done

