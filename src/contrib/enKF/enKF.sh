#!/usr/bin/env bash
#
# Copyright (C) 2017-2026, Marco Bajo, CNR-ISMAR Venice, All rights reserved.
#
# ------------------------------------------------------------------------------
# Ensemble Kalman Filter (EnKF) for SHYFEM
# ------------------------------------------------------------------------------
#
# See the README file for an help
#
#


# --- PATH & ENVIRONMENT SETUP ---
SCRIPT=$(realpath $0)
SCRIPTPATH=$(dirname $SCRIPT)
SRCDIR=$SCRIPTPATH/../..	# src directory
SIMDIR=$(pwd)		# current dir

#----------------------------------------------------------
Usage() {
    echo "Usage: enKF.sh [method] [localisation] [n] [out]"
    echo ""
    echo "Arguments:"
    echo "  method        : Analysis algorithm (11|12|13|21|22|23)"
    echo "  localisation  : Spatial localization (0: Off, 1: On)"
    echo "  n             : Threads/cores"
    echo "  out           : Output (0: mean/std only, 1: save all restarts)"
    exit 1
}

#----------------------------------------------------------
Check_file() {
    if [ ! -s "$1" ]; then
        echo "[ERROR] File missing or zero size: $1"
        exit 1
    fi
}

#----------------------------------------------------------
Check_files() {
    echo "[INFO] Validating executables..."
    command -v parallel > /dev/null 2>&1 || { echo "[ERROR] GNU Parallel not found."; exit 1; }
    [ ! -s "$SRCDIR/shyfem/shyfem" ] && echo "[ERROR] SHYFEM binary missing." && exit 1
    [ ! -s "$SRCDIR/shyfem/rstinf" ] && echo "[ERROR] SHYFEM binary missing." && exit 1
    [ ! -s "$SRCDIR/contrib/enKF/main" ] && echo "[ERROR] EnKF main binary missing." && exit 1

    echo "[INFO] Validating input lists..."
    for f in ens_list.txt obs_list.txt antime_list.txt; do Check_file "$f"; done

    echo "[INFO] Validating bas file..."
    bas_file=$(ls *.bas | head -1 2>/dev/null)
    [[ -z "$bas_file" ]] && echo "[ERROR] bas file missing." && exit 1


    rstfile_init=$(cat ens_list.txt | awk '{print $2; exit}')
    Check_file "$rstfile_init"
    echo "[INFO] Grabbing dimensions from the restart file..."
    nnkn=$($SRCDIR/shyfem/rstinf $rstfile_init | awk '/nkn +nel +nlv/ {getline; print $3; exit}')
    nnel=$($SRCDIR/shyfem/rstinf $rstfile_init | awk '/nkn +nel +nlv/ {getline; print $4; exit}')
    nnlv=$($SRCDIR/shyfem/rstinf $rstfile_init | awk '/nkn +nel +nlv/ {getline; print $5; exit}')
    echo "[INFO] Dimensions: $nnkn, $nnel, $nnlv"
}

#----------------------------------------------------------
Read_ens_list() {
    echo "[INFO] Initializing ensemble members..."
    rm -f an00001_en*b.rst
    rst1=$(head -1 ens_list.txt | awk '{print $2}')
    rst2=$(head -2 ens_list.txt | tail -1 | awk '{print $2}')
    [[ "$rst1" = "$rst2" ]] && echo "Warning!!! Perturbe the RST file!"

    nrow=0
    while read -r skelf rstf || [ -n "$skelf" ]; do
        [ -z "$skelf" ] && continue
        Check_file "$skelf"
        Check_file "$rstf"

        skel_file[$nrow]="$skelf"

        nel=$(printf "%05d" "$nrow")
        ln -fs "$rstf" "an00001_en${nel}b.rst"

        nrow=$((nrow + 1))
    done < ens_list.txt
    nrens=$nrow
}

#----------------------------------------------------------
Read_antime_list() {
    nrow=0; nran=0
    while read -r line || [ -n "$line" ]; do
        nrow=$((nrow + 1))
        if [ $((nrow % 2)) -eq 1 ]; then
            nran=$((nran+1))
            timeo[$nran]=$(echo "$line" | awk '{print $1}')
            nfile[$nran]=$(echo "$line" | awk '{print $2}')
        else
            isfile[$nran]="$line"
        fi
    done < antime_list.txt
}

#----------------------------------------------------------
SkelStr() {
    sed -e "s|NAMESIM|$1|g" -e "s|ITANF|$2|g" \
        -e "s|ITEND|$3|g" -e "s|RESTRT|$4|g" \
        -e "s|IDTRST|-1|g" "$5" > "$6"
}

#----------------------------------------------------------
Write_obs_file() {
    local na=$1
    IFS=' ' read -r -a nisfile <<< "${isfile[$na]}"
    rm -f obs_list_tmp.txt
    local k=0
    while read -r line || [ -n "$line" ]; do
        if [ "${nisfile[$k]}" = "1" ]; then echo "$line" >> obs_list_tmp.txt; fi
        k=$((k+1))
    done < obs_list.txt
}

#----------------------------------------------------------
Write_info_file() {
    local na=$1
    {
        echo "$nnkn"; echo "$nnel"; echo "$nnlv"; 
        echo "$nrens"; echo "$na"; echo "$bas_file"
        echo "${timeo[$na]}"; echo "obs_list_tmp.txt"
        echo "$rmode"; echo "$islocal"
    } > analysis.info
}

#----------------------------------------------------------
Run_ensemble_analysis() {
    local na=$1
    local nanl=$(printf "%05d" "$na")
    echo "[ANALYSIS] Executing Fortran EnKF..."
    "$SRCDIR/contrib/enKF/main"
    [ $? -ne 0 ] && echo "[ERROR] EnKF Core failed." && exit 1

    # Final check for analysis restarts
    for (( ne = 0; ne < nrens; ne++ )); do
        nensl=$(printf "%05d" "$ne")
	filename="an${nanl}_en${nensl}a.rst"
	Check_file $filename
    done
}

# -------------------------------------------------------------------
# ------------------------------ MAIN -------------------------------
# -------------------------------------------------------------------

[ $# -ne 4 ] && Usage
rmode=$1; islocal=$2; nthreads=$3; out_verb=$4

# Checking the executable programs
Check_files

# Reading skel file list
Read_ens_list

# Reading obs file list
Read_antime_list

# Assimilation cycle for every analysis time step
echo "running Assimilation cycle..."

rm -f X5*.* X3*.* backKF_*.rst analKF_*.rst
export OMP_NUM_THREADS=$nthreads
for (( na = 1; na <= nran; na++ )); do
   echo -e "\n--- Assimilation cycle STEP $na OF $nran ---"

   Write_obs_file "$na"

   Write_info_file "$na"

   # 1. ANALYSIS
   Run_ensemble_analysis "$na"

   # 2. FORECAST (only if not the last step)
   if [ "$na" -ne "$nran" ]; then
      echo "[FORECAST] Advancing ensemble... $na/$nran"
      str_list=""
      for (( ne = 0; ne < nrens; ne++ )); do
         nel=$(printf "%05d" "$ne"); nal=$(printf "%05d" "$na")
         naa=$((na + 1)); naal=$(printf "%05d" "$naa")
         name_sim="an${naal}_en${nel}b"
         strname="${name_sim}.str"
         SkelStr "$name_sim" "${timeo[$na]}" "${timeo[$naa]}" "an${nal}_en${nel}a.rst" "${skel_file[$ne]}" "$strname"
         str_list="$str_list $strname"
      done

      export OMP_NUM_THREADS=1
      # This saves the error files, if they are not present in the same moment
      parallel --jobs "$nthreads" "
        $SRCDIR/shyfem/shyfem {} > {.}.log 2>&1
        if [ -f fort.999 ]; then
          mv fort.999 fort.999_{.}
          echo 'Process {} failed: fort.999 saved as fort.999_{.}'
	  exit
        fi
        " ::: $str_list
      #parallel --jobs "$nthreads" "$SRCDIR/shyfem/shyfem {} > {.}.log 2>&1" ::: $str_list
      export OMP_NUM_THREADS=$nthreads
   fi

   # merge the rst files
   nanl=$(printf "%05d" $na)
   for (( ne = 0; ne < $nrens; ne++ )); do
        nel=$(printf "%05d" $ne)
        filename1="an${nanl}_en${nel}b.rst"
        filename2="an${nanl}_en${nel}a.rst"
        Check_file $filename1
        Check_file $filename2
	# If not verbose, saves only the last rst for each member
	if [ "$out_verb" -eq "1" ]; then
		# This is not used:
		#[[ "$na" -gt "1" ]] && cat $filename1 >> backKF_en$nel.rst
		cat $filename2 >> analKF_en$nel.rst
	else
	        [[ "$na" -eq "$nran" ]] && mv -f $filename2 analKF_en$nel.rst
	fi
	rm -f $filename1 $filename2
   done
   filename1="an${nanl}_mean_a.rst"
   filename2="an${nanl}_std_a.rst"
   Check_file $filename1
   Check_file $filename2
   cat $filename1 >> analKF_mean.rst
   cat $filename2 >> analKF_std.rst
   rm -f $filename1 $filename2
   rm -f an*_en*b.inf an*_en*.log an*_en*b.str 
done
echo -e "\n[SUCCESS] All files saved in the current directory."
