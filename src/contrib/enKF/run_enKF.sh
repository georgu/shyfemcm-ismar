#!/bin/bash
# This script runs the EnOI or the EnKF or the EnKS

########## INIT PARS #############
# Definition of paths
FEMDIR="`pwd`/.."
bindir=$FEMDIR/enKF

# Number of processes to run in parallel
nprocesses=8

# Name of the initial str file
strname_orig='nador.str'


#######################################################################
#######################################################################
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
rstname=$4
forc_type=$5
forcfile=$6

if [ $strname = 'ensstr_sim0.str' ]; then
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

else
   it_check=$(echo "($itend <= $itanf)" | bc)
   if [[ $it_check != 0 ]]; then
      echo "Bad time values: $itanf $itend"
      exit 1
   fi

   itrst=$itanf
   $FEMDIR/fembin/strparse.pl -value=idtrst -replace=-1 $assdir/$strname_orig
   $FEMDIR/fembin/strparse.pl -value=itrst -replace=$itrst replace.str
   #$FEMDIR/fembin/strparse.pl -value=restrt -replace=\'$rstname\' replace.str #Put in para
   $FEMDIR/fembin/strparse.pl -value=itanf -replace=$itanf replace.str
   $FEMDIR/fembin/strparse.pl -value=itend -replace=$itend replace.str
   #$FEMDIR/fembin/strparse.pl -value=$forc_type -replace=$forcfile replace.str #Put in para
   #### TMP
   namesim=$(basename $strname .str)
   cat replace.str | sed -e "s/Nador_sim1/$namesim/" | \
   sed -e "/^\$name/a restrt='$rstname'" | \
   sed -e "s/$forc_type.\+/$forc_type = \'input\/$forcfile\'/" > replace1.str
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
rstname=$3
forc_type=$4
forc_basename=$5

itrst=$itanf

for ((ne = 1; ne <= $nrens; ne++)); do
   echo "Initial ensemble run n. $ne of $nrens"
   nelab=$(printf "%04d" $ne)

   forcfile="${forc_basename}${nelab}.dat"
   strname="ensstr_ob0001_ens${nelab}.str"

   make_str $itanf $itend $strname $rstname "boundn" "$forcfile"

   cd $assdir
   ./shyfem $strname > ${ne}.log &
   cd -

   # set in parallel
   nemod=$((ne%$nprocesses))
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
for (( ne = 1; ne <= $nrens; ne++ )); do
    nelab=$(printf "%04d" $ne)
    rstname="ensstr_ob${nolab}_ens${nelab}.rst"
    if [ ! -s $assdir/$rstname ]; then
       echo "Bad rst file: $assdir/$rstname" 
       exit 1
    fi
    echo "$rstname" >> $assdir/rstfile_list.txt
done

cd $bindir
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

# read a file with the assimilation parameters
read_asspar

cd $assdir
rm -f shyfem; ln -s $FEMDIR/fem3d/shyfem


########## PREPARE OBS #############
# Make something to have the initial observation files with a right
# format and the obs_list.txt file with [ obs_type name_file ]
# TODO : something more general
# Creates a obs list file
echo "times_lev  obs_nador.dat" > obs_list.txt	#TMP

# Create just one binary file with all the observations
cd $bindir
make make_obs_file
./make_obs_file
cd -

# Read the times of observations
read_obs_times

isass=0

if [ $isass != 1 ]; then # begin isass ****

########## SIM BEFORE THE ENSEMBLE CREATION #############
#
# First normal spin-up simulation

itanf=-999
it_ens_spinup=86400
itend=$(echo "${otime[1]} - $it_ens_spinup" | bc -l)
strfile="ensstr_sim0.str"
#
make_str $itanf $itend $strfile
#
cd $assdir
./shyfem ensstr_sim0.str
cd -



########## CREATION OF THE INITIAL ENS #############
## Creation of an initial ensemble of states
# Introduce a perturbation on the forcing and
# save the final ensemble at the first obs, as
# initial ensemble.

itanf=$itend			#start
itend=${otime[1]}		#end
rstname="ensstr_sim0.rst"
forc_type="boundn"		#perturbed forcing type
forc_basename="astro_nador2005_last"	#perturbed forcing basename

make_init_ens_forcing $itanf $itend $rstname $forc_type $forc_basename


## Improved sampling (optional but recommended)
# Create a bigger initial ensemble and then find the
# main orthogonal direction with a SVD (Evensen, 2004).
# TODO: see m_sample1D.F90


fi # end isass ****

########## ASSIMILATION #############
echo "*** Begin of Ensemble Kalman Filter analysis ***"

nobs=${#otime[@]}	# Array length
for (( no = 1; no <= $nobs; no++ )); do

    nolab=$(printf "%04d" $no)
    nolab_new=$(printf "%04d" $((no + 1)))
    ot=${otime[$no]}

    echo "*** Analysis step: $no of $nobs"
    
    # make analisis and create new restart files
    make_analysis $nolab $ot $nrens

    # Ensemble simulations
    itanf=$ot
    if [ $no = $nobs ]; then	# last step
       itend=$tt_end
    else
       itend=${otime[$(($no + 1))]}
    fi
    echo "Initial time: $itanf"
    echo "Final time: $itend"

    for (( ne = 1; ne <= $nrens; ne++ )); do
        nelab=$(printf "%04d" $ne)
        rstname="an_ensstr_ob${nolab}_ens${nelab}.rst"
        strname="ensstr_ob${nolab_new}_ens${nelab}.str"
        make_str $itanf $itend $strname $rstname "boundn" "astro_nador2005_last.dat"

        cd $assdir
        ./shyfem $strname > ${ne}.log &
        cd -
        # set in parallel
        nemod=$((ne%$nprocesses))
        [[ $nemod = 0 ]] && wait
    done
    wait

done
echo "*** End of Ensemble Kalman Filter analysis ***"
exit

########## LAST SIM #############
# load all the restart files and average

# run the last simulation after the assimilation

### END ###
