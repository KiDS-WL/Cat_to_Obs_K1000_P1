#!/bin/bash -u
#
# ----------------------------------------------------------------
# File Name:   cuillinrun_calc2pt_mocks.sh
# Author:      Catherine Heymans (cech@roe.ac.uk)
# Description: Queues a set of observables to be measured
#              This queue submission is to analyse mocks
# ----------------------------------------------------------------

#We want to keep logs of everything so we can easily
#track back to find out why something has failed

USER=`whoami`
LOGDIR=/home/$USER/KIDS_LOGS/2ptStats/
mkdir -p $LOGDIR

# The do_all script defaults can be changed as follows
#  -d /path/to/catalogues"
#  -o /path/to/results"
#  -p patch name N or S"
#  -m list of modes"
#  -v lensfit version"
#  -n ntomo number of tomographic source bins, followed by bin edges z_B(ntomo+1)"
#     default: TOMOINFO="6 0.1 0.3 0.5 0.7 0.9 1.2 2.0"
#  -t nbins number of theta bins, theta_min, theta_max"
#     default: BININFO="9 0.5 300"
#  -i cross correlate bins i with j"
#  -j cross correlate bins i with j"
#  -c c-corr on? true/false"
#     default: true
#  -l linear not log bins? true/false"
#     default: false
#  -b which blind?"
#  -u user defined catalogue? what is the headname?"

#The modes that you can run are
# CREATETOMO\": Cut KiDS catalogues into tomographic bins and calculate and subtract c-term"
# XI: calculate xi+/- for tomo bin combination i j"
# COSEBIS: calculate En/Bn for tomo bin combination i j "
# Pkk: calculate cosmic shear Band powers for tomo bin combination i j 
# GAMMAT: calculate gamma_t and gamma_x for cross bin combination i j"
# Pgk: calculate GGL Band powers to cross bin combination i j
# COMBINE: combine the results from N and S for cross bin combination i j"

#====================================================
# MOCK SPECIFIC INFO
#====================================================
# In the data we split the catalogue into the user defined tomographic bins
# For the mocks, catalogues are created for each tomographic bin - i.e there is
# no redshift with which to split the catalogues using the CREATETOMO mode
# We want a single script for data and mocks and so we rely on the file name of the 
# mocks to have the same format as the output from CREATETOMO which is
# someheadname_6Z_1.fits  for 6 bins in total - this is tomographic bin 1
# someheadname_6Z_2.fits  for 6 bins in total - this is for tomographic bin 2 etc
# "someheadname" is defined with the -u option

# Linc's mocks currently do not have this naming system so building in a hack in this
# run script to make symbolic links.   We need to change Linc's naming system
#The location of the mocks
LincDIR=/disk05/calin/91_Data/mockFootprint/upupa/MFP_galCat
#Where I want the output to appear
ODIR=/home/cech/TMPLDIR/

#Lets submit the XI-calculation to different nodes	  
# You can use this switch to switch on/off parts of the script
#: <<'=runmelater'
for los in {0..0} #999}
do

  for type in {0..0} #4}
  do

    echo "I will submit los and type" $los $type

    # hack here to rename Linc's catalogues
    # the 'someheadname' for these catalogues is mockcat: 
    mockcat=galCat_run${los}_type${type}
    # this is a single redshift bin analysis so add _1Z_1 to the name to make this
    # consistent with the KiDS-path in the script and place it in the TOMOCATS bin of 
    # the outdir $ODIR
    ln -s $LincDIR/$mockcat.fits  $ODIR/TOMOCATS/${mockcat}_1Z_1.fits

    #here I don't really care which worker it goes to as Treecorr will use all the
    #CPU for whichever worker I sent the script to (use --exclusive)
    #I want a datadisk as I might need to copy the catalogue to the local worker 
    #lets see what the I/O is like

    # This example is running XI for 300 theta bins between 0.1 < theta < 300 (-t "300 0.1 300")
    # for a tomographic analysis with only 1 bin spanning 0.0<ZB<3.0 (-n "1 0.0 3.0")
    # for tomographic bin combination 1:1 (-i 1 -j 1)

    # If you wanted to run XI for 10 theta bins between 0.5 < theta < 300 (-t "10 0.5 300")
    # for a tomographic analysis with 5 bins spanning 0.1<ZB<1.1 (-n "6 0.1 0.3 0.5 0.7 0.9 1.1")
    # for tomographic bin combination 3:5 (-i 3 -j 5)

    sbatch --job-name=upupa_${los}_${type} \
	        --output=$LOGDIR/XI_${los}_${type}.log \
	        --error=$LOGDIR/XI_${los}_${type}.error \
	        --time=0:30:00 \
	        --exclusive \
	        --constraint="datadisk" \
	        --tasks-per-node=1 \
	        --mem=0G \
	        doall_calc2pt.sh -m XI -u $mockcat -o $ODIR -t "300 0.1 300" -n "1 0.0 3.0" -i 1 -j 1
  done
done
#=runmelater

#Once all the 2pt results are computed you might want to calculate Pkk
: <<'=runmelater'
# Do you want to calculate Pkk?  Safe to do this on the head node
for los in {1000..1200}
do
  for type in {0..1} 
  do
    mockcat=galCat_run${los}_type${type}.fits
    ./doall_calc2pt.sh -m Pkk -u $mockcat -d $LDIR -o $ODIR -t "300 0.06 300" -n "1 0.0 3.0" -i 1 -j 1
  done
done   
=runmelater