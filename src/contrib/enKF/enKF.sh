#!/bin/bash
#
# Copyright (C) 2017-2026, Marco Bajo, CNR-ISMAR Venice, All rights reserved.
#
# ------------------------------------------------------------------------------
# Ensemble Kalman Filter (EnKF) for SHYFEM
# ------------------------------------------------------------------------------
#

SCRIPT=$(realpath $0)
SCRIPTPATH=$(dirname $SCRIPT)
SRCDIR=$SCRIPTPATH/../..	# Source code base directory
SIMDIR=$(pwd)		# Current run execution workspace directory

# Parallel simulations scaling parameters
# MPI MODE: Allocates discrete MPI processes per member execution
CORES_PER_MEMBER=1
# OPENMP MODE: Launches concurrent members utilizing localized multi-threading
THREADS_PER_MEMBER=1

#----------------------------------------------------------
# Usage Information Function
#----------------------------------------------------------
Usage() {
    echo "Usage: enKF.sh [method] [localisation] [n] [out] [parallel_mode]"
    echo ""
    echo "Arguments:"
    echo "  method        : Analysis algorithm (11|12|13|21|22|23)"
    echo "  localisation  : Spatial localization (0: Off, 1: On)"
    echo "  n             : Total threads/cores available on the machine"
    echo "  out           : Output verbosity (0: mean/std only, 1: save all restarts)"
    echo "  parallel_mode : Forecast execution engine (mpi|omp)"
    exit 1
}

#----------------------------------------------------------
# Single File Integrity Validation
#----------------------------------------------------------
Check_file() {
    if [ ! -s "$1" ]; then
        echo "[ERROR] File missing or zero size: $1"
        exit 1
    fi
}

#----------------------------------------------------------
# Executables and Inputs Validation
#----------------------------------------------------------
Check_files() {
    echo "[INFO] Validating executables..."
    command -v parallel > /dev/null 2>&1 || { echo "[ERROR] GNU Parallel not found."; exit 1; }
    [ ! -s "$SRCDIR/shyfem/shyfem" ] && echo "[ERROR] SHYFEM binary missing." && exit 1
    [ ! -s "$SRCDIR/shyfem/rstinf" ] && echo "[ERROR] SHYFEM binary missing." && exit 1
    [ ! -s "$SRCDIR/contrib/enKF/main" ] && echo "[ERROR] EnKF main binary missing." && exit 1

    echo "[INFO] Validating input lists..."
    for f in ens_list.txt obs_list.txt antime_list.txt; do Check_file "$f"; done

    echo "[INFO] Validating grid configuration (.bas) file..."
    bas_file=$(ls *.bas | head -1 2>/dev/null)
    [[ -z "$bas_file" ]] && echo "[ERROR] bas file missing." && exit 1

    # Grab domain spatial dimensions from the first available restart file
    rstfile_init=$(cat ens_list.txt | awk '{print $2; exit}')
    Check_file "$rstfile_init"
    echo "[INFO] Grabbing dimensions from the restart file..."
    nnkn=$($SRCDIR/shyfem/rstinf $rstfile_init | awk '/nkn +nel +nlv/ {getline; print $3; exit}')
    nnel=$($SRCDIR/shyfem/rstinf $rstfile_init | awk '/nkn +nel +nlv/ {getline; print $4; exit}')
    nnlv=$($SRCDIR/shyfem/rstinf $rstfile_init | awk '/nkn +nel +nlv/ {getline; print $5; exit}')
    echo "[INFO] Dimensions parsed: Nodes=$nnkn, Elements=$nnel, Levels=$nnlv"
}

#----------------------------------------------------------
# Ensemble Members Parsing and Symlinking
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
# Analysis Timesteps List Parsing
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
# Simulation Skeleton Input Customization (sed replacement)
#----------------------------------------------------------
SkelStr() {
    local pefile="$6"
    local sed_args=()
    local col1 col2 dummy

    if [ -s "$pefile" ]; then
        while read -r col1 col2 dummy || [ -n "$col1" ]; do
            [[ -z "$col1" || -z "$col2" || "$col1" == \#* ]] && continue
            sed_args+=("-e" "s|${col1}|${col2}|g")
        done < "$pefile"
    fi

    sed "${sed_args[@]}" \
        -e "s|NAMESIM|$1|g" \
        -e "s|ITANF|$2|g" \
        -e "s|ITEND|$3|g" \
        -e "s|RESTRT|$4|g" \
        -e "s|IDTRST|-1|g" \
        "$5" > "$7"
}

#----------------------------------------------------------
# Filter Active Observations for the Current Timestep
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
# Generate Metadata Runtime Configuration File for EnKF Core
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
# Fortran EnKF Core Analysis Execution Block
#----------------------------------------------------------
Run_ensemble_analysis() {
    local na=$1
    local nanl=$(printf "%05d" "$na")
    echo "[ANALYSIS] Executing Fortran EnKF..."
    "$SRCDIR/contrib/enKF/main"
    [ $? -ne 0 ] && echo "[ERROR] EnKF Core failed." && exit 1

    for (( ne = 0; ne < nrens; ne++ )); do
        nensl=$(printf "%05d" "$ne")
        filename="an${nanl}_en${nensl}a.rst"
        Check_file $filename
    done
}

#==========================================================
# MAIN EXECUTION CORE
#==========================================================
# Ensure all 5 required input parameters are parsed
[ $# -ne 5 ] && Usage
rmode=$1; islocal=$2; nthreads=$3; out_verb=$4; parallel_mode=$5

# --- DYNAMIC HARDWARE OVERLOAD PROTECTION ---
# Automatically detect total physical CPU cores available on the system
SYSTEM_CORES=$(lscpu -p | grep -v '^#' | sort -u -t, -k 2,4 | wc -l)
[ -z "$SYSTEM_CORES" ] || [ "$SYSTEM_CORES" -le 0 ] && SYSTEM_CORES=$(nproc)

# Reserve 2 cores for OS and asynchronous I/O operations to guarantee host stability
SAFE_CORES_LIMIT=$(( SYSTEM_CORES - 2 ))
[ $SAFE_CORES_LIMIT -lt 2 ] && SAFE_CORES_LIMIT=2

# Throttle execution threads down to the safe threshold if over-allocated by user
if [ "$nthreads" -gt "$SAFE_CORES_LIMIT" ]; then
    echo "[WARNING] Requested $nthreads threads, but the safe limit for this machine is $SAFE_CORES_LIMIT cores."
    echo "[WARNING] Automatically adjusting 'nthreads' to $SAFE_CORES_LIMIT to prevent MPI/RAM starvation."
    nthreads=$SAFE_CORES_LIMIT
fi

# Initialize file structures, list dimensions, and global arrays mapping
Check_files
Read_ens_list
Read_antime_list

# --- AUTOMATIC RESTRT INTEGRITY & SYNCHRONIZATION (SPIN-UP) ---
echo "[INFO] Verifying if first observation record (${timeo[1]}) matches or exists in the initial restart..."

# Query rstinf output via grep to check if the target assimilation timestamp is already present in the file history
date_exists=$($SRCDIR/shyfem/rstinf $rstfile_init | grep -F "${timeo[1]}")

if [ -n "$date_exists" ]; then
    echo "[SUCCESS] First observation record (${timeo[1]}) is natively available in the restart history."
    echo "[SUCCESS] SHYFEM will parse the index correctly. Skipping synchronization spin-up."
else
    echo "[WARNING] First observation record (${timeo[1]}) NOT found in the restart file."
    # Extract the final calculated date record string from the initial restart file to use as the spin-up baseline
    rst_date_final=$($SRCDIR/shyfem/rstinf $rstfile_init | awk '/Final time in file/ {print $6}')
    echo "[SYNCHRONIZATION] Advancing all $nrens members from $rst_date_final to ${timeo[1]}..."

    # 1. Generate synchronization skeleton files (.str) for all ensemble members
    for (( ne = 0; ne < nrens; ne++ )); do
        nel=$(printf "%05d" "$ne")
        name_sim="sync_en${nel}b"
        pefile="pe_parameters_init_en${nel}.dat" 
        strname="${name_sim}.str"
        
        # Invoke SkelStr setting idtrst=-1. SHYFEM will output a clean, single-record restart at ITEND (timeo[1])
        SkelStr "$name_sim" "$rst_date_final" "${timeo[1]}" "an00001_en${nel}b.rst" "${skel_file[$ne]}" "$pefile" "$strname"
    done

    # 2. Execute the parallel spin-up computational block (MPI vs OpenMP)
    if [ "$parallel_mode" = "mpi" ]; then
        export OMP_NUM_THREADS=1
        JOBS_CONCURRENT=$(( nthreads / CORES_PER_MEMBER ))
        [ $JOBS_CONCURRENT -lt 1 ] && JOBS_CONCURRENT=1
        runpr=$([ $CORES_PER_MEMBER -eq "1" ] && echo "" || echo "mpirun -np $CORES_PER_MEMBER")

        parallel --halt now,fail=1 --jobs "$JOBS_CONCURRENT" "
            $runpr $SRCDIR/shyfem/shyfem {} > {.}.log 2>&1 || true
            if ! grep -q '100.000 %' {.}.log; then
                [ -f fort.999 ] && mv fort.999 fort.999_{.}
                exit 1
            else
                [ -f fort.999 ] && rm -f fort.999; exit 0
            fi
        " ::: sync_en*b.str
    else
        export OMP_NUM_THREADS=$THREADS_PER_MEMBER
        JOBS_CONCURRENT=$(( nthreads / THREADS_PER_MEMBER ))
        [ $JOBS_CONCURRENT -lt 1 ] && JOBS_CONCURRENT=1

        parallel --halt now,fail=1 --jobs "$JOBS_CONCURRENT" "
            $SRCDIR/shyfem/shyfem {} > {.}.log 2>&1 || true
            if ! grep -q '100.000 %' {.}.log; then
                [ -f fort.999 ] && mv fort.999 fort.999_{.}
                exit 1
            else
                [ -f fort.999 ] && rm -f fort.999; exit 0
            fi
        " ::: sync_en*b.str
    fi

    # Evaluate the global exit status of the synchronization execution
    if [ $? -ne 0 ]; then
        echo "[ERROR] Spin-up alignment failed. Check sync_en*.log execution outputs."
        exit 1
    fi

    # 3. Splice newly generated single-record restarts into the execution core workspace
    echo "[SYNCHRONIZATION] Splicing single-record synchronized restarts into execution core..."
    for (( ne = 0; ne < nrens; ne++ )); do
        nel=$(printf "%05d" "$ne")
        Check_file "sync_en${nel}b.rst"

        # Replace the initial symbolic link with the actual advanced single-record restart file
        rm -f "an00001_en${nel}b.rst"
        mv -f "sync_en${nel}b.rst" "an00001_en${nel}b.rst"
        rm -f sync_en${nel}b.log sync_en${nel}b.str
    done
    echo "[SUCCESS] Spin-up completed. Workspace successfully aligned to ${timeo[1]}."
fi

echo "Starting Assimilation Cycle..."

# Purge leftover temporary structures from previous runs
rm -f X5*.* X3*.* backKF_*.rst analKF_*.rst

# Maximize multi-threading allocation for the initialization phase of the EnKF compiled core
export OMP_NUM_THREADS=$nthreads

# --- MAIN ASSIMILATION LOOP ---
for (( na = 1; na <= nran; na++ )); do
    echo -e "\n--- Assimilation Cycle STEP $na OF $nran ---"
    Write_obs_file "$na"
    Write_info_file "$na"

    # 1. ANALYSIS STEP (Symmetric OpenMP multiprocessing within the Fortran binary core)
    Run_ensemble_analysis "$na"

    # 2. FORECAST STEP (Evaluated exclusively if a subsequent timestep is pending)
    if [ "$na" -ne "$nran" ]; then
        echo "[FORECAST] Advancing ensemble members... $na/$nran"

        # Generate custom simulation string configurations (.str) for each ensemble member
        for (( ne = 0; ne < nrens; ne++ )); do
            nel=$(printf "%05d" "$ne"); nal=$(printf "%05d" "$na")
            naa=$((na + 1)); naal=$(printf "%05d" "$naa")
            name_sim="an${naal}_en${nel}b"
            pefile="pe_parameters_an${nal}_en${nel}.dat"
            strname="${name_sim}.str"
            SkelStr "$name_sim" "${timeo[$na]}" "${timeo[$naa]}" "an${nal}_en${nel}a.rst" "${skel_file[$ne]}" "$pefile" "$strname"
        done

        # --- DYNAMIC FORECAST PARALLELIZATION PARSING ---
        if [ "$parallel_mode" = "mpi" ]; then
            # MPI Mode execution parameters mapping
            export OMP_NUM_THREADS=1
            JOBS_CONCURRENT=$(( nthreads / CORES_PER_MEMBER ))
            [ $JOBS_CONCURRENT -lt 1 ] && JOBS_CONCURRENT=1
            if [ $CORES_PER_MEMBER -eq "1" ]; then runpr='' ; else runpr="mpirun -np $CORES_PER_MEMBER" ; fi

            echo "[INFO] [MPI MODE] Running $JOBS_CONCURRENT concurrent members via $runpr..."
            parallel --halt now,fail=1 --jobs "$JOBS_CONCURRENT" "
                $runpr $SRCDIR/shyfem/shyfem {} > {.}.log 2>&1 || true
                if ! grep -q '100.000 %' {.}.log; then
                    [ -f fort.999 ] && mv fort.999 fort.999_{.}
                    echo 'Process {} failed (100% not reached).'
                    exit 1
                else
                    [ -f fort.999 ] && rm -f fort.999
                    exit 0
                fi
            " ::: an${naal}_en*b.str

            if [ $? -ne 0 ]; then
                echo "[ERROR] Forecast step failed in MPI execution. Check the remaining .log files."
                exit 1
            fi
        else
            # OpenMP Mode execution parameters mapping
            export OMP_NUM_THREADS=$THREADS_PER_MEMBER
            JOBS_CONCURRENT=$(( nthreads / THREADS_PER_MEMBER ))
            [ $JOBS_CONCURRENT -lt 1 ] && JOBS_CONCURRENT=1

            echo "[INFO] [OMP MODE] Running $JOBS_CONCURRENT concurrent members, each restricted to $THREADS_PER_MEMBER OpenMP threads..."
            parallel --halt now,fail=1 --jobs "$JOBS_CONCURRENT" "
                $SRCDIR/shyfem/shyfem {} > {.}.log 2>&1 || true
                if ! grep -q '100.000 %' {.}.log; then
                    [ -f fort.999 ] && mv fort.999 fort.999_{.}
                    echo 'Process {} failed (100% not reached).'
                    exit 1
                else
                    [ -f fort.999 ] && rm -f fort.999
                    exit 0
                fi
            " ::: an${naal}_en*b.str

            if [ $? -ne 0 ]; then
                echo "[ERROR] Forecast step failed in OMP execution. Check the remaining .log files."
                exit 1
            fi
        fi

        # Restore maximum available threads for the subsequent EnKF Core analysis execution
        export OMP_NUM_THREADS=$nthreads
    fi

    # --- CONSOLIDATE AND MERGE RESTART OUTPUT ---
    nanl=$(printf "%05d" $na)
    for (( ne = 0; ne < $nrens; ne++ )); do
        nel=$(printf "%05d" $ne)
        filename1="an${nanl}_en${nel}b.rst"
        filename2="an${nanl}_en${nel}a.rst"
        Check_file $filename1
        Check_file $filename2

        # Output verbosity check: append complete restart history or isolate the last time-slice
        if [ "$out_verb" -eq "1" ]; then
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

    # SELECTIVE CLEANUP: Clean up configuration logs and structures if current step succeeded
    if [ "$na" -ne "$nran" ]; then
        rm -f an${naal}_en*b.log an${naal}_en*b.str an${naal}_en*b.inf
    fi
done

echo -e "\n[SUCCESS] Assimilation cycle complete. All final assets consolidated in the current workspace."
