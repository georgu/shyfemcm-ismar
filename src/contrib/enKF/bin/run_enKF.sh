#!/bin/bash
# This script runs the EnOI or the EnKF or the EnKS

#####
# Reads a file with the input parameters
read_asspar(){
assdir=$(sed -n "1p" "asspar.txt")	# working dir
nrens=$(sed -n "2p" "asspar.txt")	# nr ens members
adate=$(sed -n "3p" "asspar.txt")	# date at time 0: yyyymmdd
atime=$(sed -n "4p" "asspar.txt")	# time
tt_anf=$(sed -n "5p" "asspar.txt")	# first time
tt_end=$(sed -n "6p" "asspar.txt")	# last time
tt_eanf=$(sed -n "7p" "asspar.txt")	# first time ens members
tt_eend=$(sed -n "8p" "asspar.txt")	# last time ens members
tt_oanf=$(sed -n "9p" "asspar.txt")	# first time obs (analysis)
tt_oend=$(sed -n "10p" "asspar.txt")	# last time obs (analysis)
}

#####
# Read the times contained in the obs_times.dat and store
# them into an array
read_obs_times(){
obsfile=$assdir/obs_times.dat
k=0
while read line
do
    k=$((k+1))
    otime[$k]=$line
done < $obsfile
}

#####
# Make the str files
make_str(){
itanf=$1
itend=$2
strname=$3

if [ $strname = 'ensstr_sim0.str' ]; then
   echo "Initial run"
   itanf0=$($FEMDIR/fembin/strparse.pl -value=itanf $assdir/$strname_orig | tail -1)
   itend0=$($FEMDIR/fembin/strparse.pl -value=itend $assdir/$strname_orig | tail -1)
   itanf=$itanf0
   it_check=$(echo "(($itend <= $itanf) || ($itend >= $itend0))" | bc)
   if [[ $it_check != 0 ]]; then
      echo "Bad time values: $itanf $itend $itend0"
      exit 1
   fi
   $FEMDIR/fembin/strparse.pl -value=idtrst -replace=-1 $assdir/$strname_orig
   $FEMDIR/fembin/strparse.pl -value=itanf -replace=$itanf replace.str
   $FEMDIR/fembin/strparse.pl -value=itend -replace=$itend replace.str
   #### TMP
   namesim=$(basename $strname .str)
   cat replace.str | sed -e "s/Nador_sim1/$namesim/" > replace1.str
   #### TMP
   mv replace1.str $assdir/$strname
   rm replace.str
   
fi
}

#####
# Run nrens simulations with different forcing to create an initial ensemble
make_init_ens_forcing(){
itanf=$1
itend=$2
itrst=$itanf
strfile=$3
nrens=$4
forc_type=$5
forc_basename=$6

namerst=$(basename $strfile .str); namerst=${namerst}.rst
for ((ne = 1; ne <= $nrens; ne++)); do
   echo "Initial ensemble run n. $ne of $nrens"
   nelab=$(printf "%04d" $ne)

   forcfile=${forc_basename}${nelab}.dat

   $FEMDIR/fembin/strparse.pl -value=idtrst -replace=-1 $assdir/$strname_orig
   $FEMDIR/fembin/strparse.pl -value=itrst -replace=$itrst replace.str
   #$FEMDIR/fembin/strparse.pl -value=restrt -replace=\'$namerst\' replace.str #Put in para
   $FEMDIR/fembin/strparse.pl -value=itanf -replace=$itanf replace.str
   $FEMDIR/fembin/strparse.pl -value=itend -replace=$itend replace.str
   #$FEMDIR/fembin/strparse.pl -value=$forc_type -replace=$forcfile replace.str #Put in para
   filename="ensstr_ob0000_ens${nelab}.str"
   title=$(basename $filename .str)
   #### TMP
   cat replace.str | sed -e "s/Nador_sim1/$title/" | \
   sed -e "/^\$name/a restrt='$namerst'" | \
   sed -e "s/$forc_type.\+/$forc_type = \'input\/$forcfile\'/" > replace1.str 
   #### TMP
   mv replace1.str $assdir/$filename

   cd $assdir
   $FEMDIR/fem3d/ht < $filename > ${ne}.log &
   cd -

   # set in parallel
   nemod=$((ne%$nprocesses)); echo "************ Number of run: $ne ****************"
   [[ $nemod = 0 ]] && wait
done
wait
}


#####
# Make the analysis at every observation time
make_analysis(){
nolab=$1
tt=$2
nrens=$3

# Make a list of rst files and check their size
rm -f $assdir/rstfile_list.txt
for ((ne = 1; ne <= $nrens; ne++)); do
    nelab=$(printf "%04d" $ne)
    rstname="ensstr_ob${nolab}_ens${nelab}.rst"
    if [ ! -s $assdir/$rstname ]; then
       echo "Bad rst file: $assdir/$rstname" 
       exit 1
    fi
    echo "$assdir/$rstname" >> $assdir/rstfile_list.txt
done

echo "$tt" > $assdir/asstime.txt
make enKF_analysis
./enKF_analysis < $assdir/asstime.txt
}

######################################################################################
#					 MAIN					     #
######################################################################################
#
# This is the assimilation scheme:
#
#        |------| |------| |------| |------| | 
#        |------| |------| |------| |------| |
#        |------| |------| |------| |------| |
#        |------| |------| |------| |------| |
#  ------|------|A|------|A|------|A|------|A|-<A>-----------
#        |------| |------| |------| |------| |
#        |------| |------| |------| |------| |
#        |------| |------| |------| |------| |
#        |------| |------| |------| |------| |
#
#  |Init.| |First| |Anal.|                 |Anal.| |             | |        |
#  |sim. |-|ens. |-|first|-     ...       -|last |-|Average state|-|Last run|
#  |     | |     | |obs  |                 |obs  | |             | |        |
#
#
######################################################################################

########## INIT PARS #############
# Definition of paths
FEMDIR=$HOME/marcob/enKF/shyfem-7_1_12
bindir=$FEMDIR/enKF/bin

# Number of processes to run in parallel
nprocesses=8

# Name of the initial str file
strname_orig='nador.str'

# read a file with the assimilation parameters
read_asspar
########## INIT PARS #############


########## PREPARE OBS #############
# Collect and sort observations (if T-S profiles or timeseries or very near measurements..)
cd $bindir
make make_obs_file
./make_obs_file
cd -

# Read the times of observations
read_obs_times
########## PREPARE OBS #############



########## SIM BEFORE THE ENSEMBLE CREATION #############
#
# First normal spin-up simulation

#itanf=-999
#itend=$(echo "${otime[1]} - $it_ens_spinup" | bc -l)
#strfile="ensstr_sim0.str"
#
#make_str $itanf $itend $strfile
#
#cd $assdir
#time $FEMDIR/fem3d/ht < ensstr_sim0.str
#cd -
########## SIM BEFORE THE ENSEMBLE CREATION #############



########## CREATION OF THE INITIAL ENS #############
## Creation of an initial ensemble of states
# Introduce a perturbation on the forcing and
# save the final ensemble at the first obs, as
# initial ensemble.

#itanf=$itend			#start
#itend=${otime[1]}		#end
#strfile="ensstr_sim0.str"	#base str file
#forc_type="boundn"		#perturbed forcing type
#forc_basename="astro_nador2005_last"	#perturbed forcing basename
#
#make_init_ens_forcing $itanf $itend $strfile $nrens $forc_type $forc_basename


## Improved sampling (optional but recommended)
# Create a bigger initial ensemble and then find the
# main orthogonal direction with a SVD (Evensen, 2004).
# TODO: see m_sample1D.F90
########## CREATION OF THE INITIAL ENS #############



########## ASSIMILATION #############
echo "***************************************************"
echo "	Starting the assimilation runs"
echo "***************************************************"
no=0	# num of obs times
# obs loop
for ot in ${otime[@]}; do

    nolab=$(printf "%04d" $no)
    echo "*** Observation time $no"
    
    ######
    # make analisis and create new restart files
    make_analysis $nolab $ot $nrens
    ######
exit

    no=$((no+1));

    # Ensemble loop
    itanf=$ot; itend=${otime[$((no+1))]}
    for ((ne = 1; ne <= $nrens; ne++)); do
        nelab=$(printf "%04d" $ne)
        filename="ensstr_ob${nolab}_ens${nelab}.str"
        #make_str $itanf $ot $filename
	# simulations &
    done

done
exit
########## ASSIMILATION #############

# load all the restart files and average

# run the last simulation after the assimilation

### END ###


