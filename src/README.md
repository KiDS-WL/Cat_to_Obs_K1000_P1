# Code to convert finely binned xi_\pm measurements into COSEBIs and Bandpower estimates
This code is called from https://github.com/KiDS-WL/Cat_to_Obs_K1000_P1/blob/master/Calc_2pt_Stats/doall_calc2pt.sh using modes COSEBIS, Pkk and Pgk.

## bandpowers (xi2bandpow.c from Benjamin Joachimi)
A single c-code with the following command line arguments
1. working directory
2. input file identifier
3. output file identifier
4. number of input angular bins
5. min input separation to use in conversion [arcmin] (xi_+ in case of ee)
6. max input separation to use in conversion [arcmin] (xi_+ in case of ee)
7. min input separation to use in conversion [arcmin] (xi_- in case of ee; otherwise unused)
8. max input separation to use in conversion [arcmin] (xi_- in case of ee; otherwise unused)
9. number of output angular bins
10. min output angular frequency
11. max output angular frequency
12. correlation type (1: ee; 2: ne; 3: gg)
13. log width of apodisation window [total width of apodised range is tmax/tmin=exp(width) in arcmin; <0 for no apodisation]

Outputs #   ell         PEE             PEEerr          PBB                PBBerr, where the errors are estimated shape-noise only and should only be used as a ballpark underestimate.

The makefile links to necessary libraries in bjutils.

## cosebis (run_measure_cosebis_cat2stats.py from Marika Asgari)
Python code with the following command line arguments
* -i Full Input file name (i.e the TreeCorr xi_pm ascii output)
* -t theta column in input file
* -p xi_plus column in input file
* -m xi_minus column in input file
* --cfoldername full name and address of the folder for En/Bn files
* -o output file name suffix
* -b binning (log or linear)
* -n number of COSEBIs modes to calculate
* -s value of thetamin in arcmins
* -l value of thetamax in arcmins
* --tfoldername name and full address of the folder for Tplus Tminus files, will make it if it does not exist
* -- tplusfile name of Tplus file, will look for it before running the code
* -- tminusfile name of Tminus file, will look for it before running the code
* -c normalisation file name and address for T_plus/minus
* -r roots file name and address for T_plus/minus

```
# example run:
# tmin=0.50
# tmax=100.00
# python run_measure_cosebis_cats2stats.py -i xi_nBins_1_Bin1_Bin1 -o nBins_1_Bin1_Bin1 --norm ./TLogsRootsAndNorms/Normalization_${tmin}-${tmax}.table -r 
# ./TLogsRootsAndNorms/Root_${tmin}-${tmax}.table -b lin --thetamin ${tmin} --thetamax ${tmax} -n 20
```

The filter functions required to convert xi_pm to COSEBIs are created the first time you run the code and stored in TLogsRootsAndNorms.   This can be slow!   After you've run this once, however, you only need to recalculate them if you decide to change the theta min/max that you wish to analyse.  
