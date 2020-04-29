#!/bin/bash -u
#
# ----------------------------------------------------------------
# File Name:   cuillinrun_V1.0.0.sh
# Author:      Catherine Heymans (cech@roe.ac.uk)
# Description: Queues a set of observables to be measured
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

#The modes that you can run are
# CREATETOMO\": Cut catalogues into tomographic bins and calculate and subtract c-term"
# XI: calculate xi+/- for tomo bin combination i j"
# COSEBIS: calculate En/Bn for tomo bin combination i j "
# Pkk: calculate cosmic shear Band powers for tomo bin combination i j 
# GT: calculate gamma_t and gamma_x for cross bin combination i j"
# Pgk: calculate GGL Band powers to cross bin combination i j
# COMBINE: combine the results from N and S for cross bin combination i j"

LENSFIT_VERSION=svn_309c_2Dbins_v2

#: <<'=runmelater'
for BLIND in A #B C
do

## Do you want to create the tomographic catalogues?  Safe to do this on the head node
# No SOM Flag Selection
#./doall_calc2pt.sh -m CREATETOMO -c true -p N -v $LENSFIT_VERSION -b $BLIND
#./doall_calc2pt.sh -m CREATETOMO -c true -p S -v $LENSFIT_VERSION -b $BLIND

# SOM Flag Selection
for FLAG_SOM in Flag_SOM_Fid #Flag_SOM_noDEEP2
do
  #./doall_calc2pt.sh -m CREATETOMO -c true -p N -v $LENSFIT_VERSION -s $FLAG_SOM -b $BLIND
  #./doall_calc2pt.sh -m CREATETOMO -c true -p S -v $LENSFIT_VERSION -s $FLAG_SOM -b $BLIND
  ./doall_calc2pt.sh -m CREATETOMO -c true -p ALL -v $LENSFIT_VERSION -s $FLAG_SOM -b $BLIND
done

done
#=runmelater

#Lets submit the XI-calculation to different nodes	    
#This is for the default 5 bins

#FLAG_SOM=false
FLAG_SOM=Flag_SOM_Fid 
#FLAG_SOM=Flag_SOM_noDEEP2

: <<'=runmelater'
for patch in N S
do
    
for ibin in {1..5}
do
	jbin=$ibin
  while [[ $jbin -le 5 ]]
  do
		echo $ibin $jbin
    echo "I will submit tomographic bins and patch" $ibin $jbin $patch

    #here I don't really care which worker it goes to as Treecorr will use all the
    #CPU for whichever worker I sent the script to (use --exclusive)
    #I want a datadisk as I might need to copy the catalogue to the local worker 
    #lets see what the I/O is like

    sbatch --job-name=XI_${patch}_${ibin}_${jbin} \
	        --output=$LOGDIR/XI_${ibin}_${jbin}_${patch}.log \
	        --error=$LOGDIR/XI_${ibin}_${jbin}_${patch}.error \
	        --time=8:40:00 \
	        --exclusive \
	        --constraint="datadisk" \
		--partition=WL \
	        --tasks-per-node=1 \
	        --mem=0G \
	        doall_calc2pt.sh -m XI -i $ibin -j $jbin -p $patch -v $LENSFIT_VERSION -s $FLAG_SOM -t "4000 0.5 300.0"
#    	        doall_calc2pt.sh -m XI -i $ibin -j $jbin -p $patch -v $LENSFIT_VERSION -s $FLAG_SOM
    ((jbin = jbin + 1))
  done
done
done

=runmelater
: <<'=runmelater'


# Do you want to combine the XI N/S results?  Safe to do this on the head node
for ibin in {1..5}
do
	jbin=$ibin
  while [[ $jbin -le 5 ]]
  do
  ./doall_calc2pt.sh -m COMBINEXI -i $ibin -j $jbin -p ALL -v $LENSFIT_VERSION -s $FLAG_SOM
  ((jbin = jbin + 1))
  done
done   
=runmelater

: <<'=runmelater'
# Do you want to calculate Pkk?  Safe to do this on the head node
for ibin in {1..5}
do
	jbin=$ibin
  while [[ $jbin -le 5 ]]
  do
  ./doall_calc2pt.sh -m Pkk -i $ibin -j $jbin -p ALL -v $LENSFIT_VERSION  -s $FLAG_SOM
  ((jbin = jbin + 1))
  done
done   

=runmelater
: <<'=runmelater'
# Do you want to calculate GGL?  do not do this on the head node
# -p N automatically connects with BOSS
# -p S automatically connects with 2dFLenS
for patch in N S
do
for ibin in {1..2}
do
	for jbin in {1..5}
  do
    echo "I will submit tomographic bins and patch" $ibin $jbin $patch
    sbatch --job-name=GT_${patch}_${ibin}_${jbin} \
	        --output=$LOGDIR/GT_${ibin}_${jbin}_${patch}.log \
	        --error=$LOGDIR/GT_${ibin}_${jbin}_${patch}.error \
	        --time=8:40:00 \
	        --exclusive \
	        --constraint="datadisk" \
	        --tasks-per-node=1 \
		--partition=WL \
	        --mem=0G \
            doall_calc2pt.sh -m GT -i $ibin -j $jbin -p $patch  -v $LENSFIT_VERSION #-s $FLAG_SOM
  done
done   
done


=runmelater
: <<'=runmelater'
# Do you want to combine the GT N/S results?  Safe to do this on the head node
for ibin in {1..2}
do
  for jbin in {1..5}
  do
    ./doall_calc2pt.sh -m COMBINEGT -i $ibin -j $jbin -p ALL -v $LENSFIT_VERSION -s $FLAG_SOM
  done
done   

=runmelater


: <<'=runmelater'
# Do you want to calculate Pgk?  Safe to do this on the head node
for ibin in {1..2}
do
    for jbin in {1..5}
    do
        ./doall_calc2pt.sh -m Pgk -i $ibin -j $jbin -p ALL -v $LENSFIT_VERSION  -s $FLAG_SOM
    done
done   
=runmelater

: <<'=runmelater'
# Do you want to calculate COSEBIS?  Safe to do this on the head node
for ibin in {1..5}
do
	jbin=$ibin
  while [[ $jbin -le 5 ]]
  do
  ./doall_calc2pt.sh -m COSEBIS -i $ibin -j $jbin -p ALL -v $LENSFIT_VERSION -s $FLAG_SOM
  ((jbin = jbin + 1))
  done
done   

=runmelater
