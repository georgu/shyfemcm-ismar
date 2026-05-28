#!/usr/bin/env bash

# --- Path Deduction ---
ENKF_DIR=$(dirname "$(readlink -f "$0")")
SHYFEM_DIR=$(cd "$ENKF_DIR/../../.." && pwd)
DET_DIR=$(pwd)
DA_DIR="$DET_DIR/sim_DA"
OBS_DIR="$DET_DIR/observations"

#===================
check_env() {
    echo ">>> Step 1: Environment Setup"
    rm -fr "$DA_DIR"
    mkdir -p "$DA_DIR"
    [ ! -f "$ENKF_DIR/enKF.sh" ] && { echo "ERROR: Main script enKF.sh not found in $ENKF_DIR."; exit 1; }
    cp -p $SHYFEM_DIR/src/external/gotm/gotmturb.nml $DA_DIR
}

#===================
ask_params() {
    echo ">>> Step 2: Ensemble Size Configuration"
    echo "This will determine the number of parallel simulations to run."
    read -p "--- [INPUT] Enter the number of ensemble members (NRENS): " NRENS
    [ $((NRENS % 2)) -eq 0 ] && NRENS=$((NRENS + 1)) && echo "--- [NOTE] NRENS adjusted to $NRENS to ensure an odd number for the EnKF algorithm."
}


#===================
perturb_forcings_interactive() {
    echo ">>> Step 3: Boundary & Forcing Perturbation"
    
    # Find the first .str file to use as a template for forcing names
    local base_str=$(ls *.str 2>/dev/null | head -n 1)
    [ -z "$base_str" ] && { echo "Error: No .str file found."; return 1; }

    # Extract unique filenames from specific sections in the .str file
    # Searches for strings between quotes ending in .dat, .fem, or .txt
    local files=$(awk '/\$end/ {in_sec=0} /\$[a-zA-Z0-9]+/ {in_sec=1} in_sec && /boundn|saltn|tempn|vel3dn|wind|rain|qflux|saltin|tempin/ {if (match($0, /[\047"][^\047"]+\.(dat|fem|txt)[\047"]/)) {print substr($0, RSTART+1, RLENGTH-2)}}' "$base_str" | sort -u)

    for f in $files; do
        # Check if the source file exists in the deterministic directory (DET_DIR)
        if [ ! -f "$DET_DIR/$f" ]; then
            echo "Warning: Forcing file $f not found in $DET_DIR. Skipping..."
            continue
        fi

        # Extract base name and extension for local operations
        local fname=$(basename "$f")
        local name="${fname%.*}"
        local ext="${fname##*.}"

        echo "------------------------------------------------"
        echo "Detected forcing file: $fname"
        read -p "--- [INPUT] Perturb this forcing to create ensemble variations? (y/n): " do_p

        if [[ "$do_p" == "y" ]]; then
            # Create a symbolic link in DA_DIR to avoid path issues during perturbation
            ln -sf "$DET_DIR/$f" "$DA_DIR/$fname"
            
            # Move into DA_DIR so perturbation scripts output results locally
            pushd "$DA_DIR" > /dev/null

            # --- CASE 1: WIND & PRESSURE (Time Series - .dat/.txt) ---
            if [[ "$fname" == *"wind"* || "$fname" == *"press"* ]] && [[ "$ext" == "dat" || "$ext" == "txt" ]]; then
                read -p "      > Wind TS Params (Format[0:u,v 1:ws,dir] Press[0:no 1:yes] STD_w STD_p Tau): " W_FMT W_PR W_STDW W_STDP W_TAU
                bash "$ENKF_DIR/perturbe_ts_wind_press.sh" "$fname" "$NRENS" "$W_FMT" "$W_PR" "$W_STDW" "$W_STDP" "$W_TAU"

            # --- CASE 2: HEAT FLUX / QFLUX (Time Series - .dat/.txt) ---
            elif [[ "$fname" == *"heat"* || "$fname" == *"qflux"* ]] && [[ "$ext" == "dat" || "$ext" == "txt" ]]; then
                read -p "      > Heat TS Params (STD_Solar STD_Temp STD_Humi STD_Cloud Tau): " H_S H_A H_H H_C H_T
                bash "$ENKF_DIR/perturbe_ts_heat.sh" "$fname" "$NRENS" "$H_S" "$H_A" "$H_H" "$H_C" "$H_T"

            # --- CASE 3: GENERIC SCALAR TIME SERIES (.dat/.txt) ---
            elif [[ "$ext" == "dat" || "$ext" == "txt" ]]; then
                read -p "      > Simple TS Params (STD MIN MAX Tau): " TS_STD TS_MIN TS_MAX TS_TAU
                "$ENKF_DIR/perturbe_ts" "$fname" "$NRENS" "$TS_STD" "$TS_TAU" "$TS_MIN" "$TS_MAX"

            # --- CASE 4: FINITE ELEMENT FILES (.fem) ---
            elif [[ "$ext" == "fem" ]]; then
                if [[ "$fname" == *"wind"* || "$fname" == *"press"* ]]; then
                    read -p "      > Wind FEM Params (Method[1:P+W, 2:W] STD Tau): " W_TYPE W_STD W_TAU
                    "$ENKF_DIR/perturbe_fem_wind_press" "$fname" "$NRENS" "$W_TYPE" "$W_STD" "$W_TAU"
                elif [[ "$fname" == *"heat"* || "$fname" == *"qflux"* ]]; then
                    read -p "      > Heat FEM Params (STD_S STD_T STD_H STD_C Tau): " H_S H_A H_H H_C H_T
                    "$ENKF_DIR/perturbe_fem_heat" "$fname" "$NRENS" "$H_S" "$H_A" "$H_H" "$H_C" "$H_T"
                else
                    read -p "      > Scalar FEM Params (STD MIN MAX Tau PType): " S_STD S_VMIN S_VMAX S_TAU S_TYP
                    "$ENKF_DIR/perturbe_fem_scalar" "$fname" "$NRENS" "$S_TYP" "$S_STD" "$S_VMIN" "$S_VMAX" "$S_TAU"
                fi
            fi

            # Clean up: remove the symbolic link, keeping the generated ensemble files
            rm "$fname"
            popd > /dev/null
        else
            # If no perturbation is needed, copy the original file as a static forcing
            echo "      - Copying original $fname as static forcing."
            cp "$DET_DIR/$f" "$DA_DIR/$fname"
        fi
    done
    rm -f $DET_DIR/newton.*
}

#===================
create_ensemble_data() {
    echo ">>> Step 4: Skeleton and Restart File Generation"
    # Identify base structure (.str) and restart (.rst) files
    local base_str=$(ls *.str 2>/dev/null | head -n 1)
    local base_rst=$(ls *.rst 2>/dev/null | head -n 1)
    local rst_name="${base_rst%.*}"

    # Extract basin name for .bas file usage
    local basin_name=$(awk '/\$title/{getline; getline; getline; print $1; exit}' "$base_str")
    local basin_f="${basin_name}.bas"
    [ -f "$DET_DIR/$basin_f" ] && cp -p "$DET_DIR/$basin_f" "$DA_DIR/"

    # --- Restart Perturbation Section ---
    if [ -f "$base_rst" ]; then
        echo "Found base restart file: $base_rst"
        read -p "--- [INPUT] Perturb initial state? (y/n): " p_rst
        if [[ "$p_rst" == "y" ]]; then
            "$SHYFEM_DIR/bin/rstinf" "$base_rst"
            read -p "    > Analysis start date (yyyy-mm-dd::HH:MM:SS): " R_DATE
            read -p "    > Vertical levels (nlv): " R_NLV
            read -p "    > Amplitude SSH/T/S (errz errt errs): " R_ERRZ R_ERRT R_ERRS
            "$ENKF_DIR/perturbe_state" "$basin_f" "$base_rst" "$R_DATE" "$R_NLV" "$NRENS" "$R_ERRZ" "$R_ERRT" "$R_ERRS"
            rm -f $DET_DIR/newton.*
        else
            read -p "--- [INPUT] Start date (yyyy-mm-dd::HH:MM:SS): " R_DATE
            for ((i=0; i<NRENS; i++)); do cp "$base_rst" "${rst_name}_$(printf "%03d" $i).rst"; done
        fi
    fi

    # --- Skeleton Generation Section ---
    echo "Generating individual member Skeletons (.skel) from $base_str..."
    for ((i=0; i<NRENS; i++)); do
        idx=$(printf "%03d" $i)
        local skel_file="member_${idx}.skel"
        local rst_file="${rst_name}_${idx}.rst"

        awk -v idx="$idx" -v da_dir="$DA_DIR" '
        BEGIN { in_para=0; in_name=0; in_title=0; title_line=0 }

        # 1. Global cleaning: remove all instances of variables we want to control or delete
        /^[ \t]*(itmrst|itrst|itanf|itend|idtrst|nocheck|restrt)[ \t]*=/ { next }

        # 2. Handle $title block (NAMESIM on the 2nd line)
        /\$title/ { in_title=1; title_line=0; print; next }
        in_title {
            title_line++
            if (title_line == 1) { print "     Sim with DA"; next }
            if (title_line == 2) { print "     NAMESIM"; next }
            if (title_line == 3) { print; in_title=0; next }
        }

        # 3. Handle $para block
        /\$para/ { 
            in_para=1; print 
            # Inject mandatory variables at the beginning of $para
            printf "        itanf  = \047ITANF\047  itend  = \047ITEND\047\n"
            printf "        itrst  = \047ITANF\047\n"
            printf "        idtrst = IDTRST\n"
            printf "        nocheck = 1\n"
            next 
        }
        
        # Inside $para, force specific timers (itmout, itmcon, itmext) to ITANF
        in_para && /^[ \t]*(itmout|itmcon|itmext)[ \t]*=/ {
            sub(/=.*/, "= \047ITANF\047"); print; next
        }

        # 4. Handle $name block
        /\$name/ { 
            in_name=1; print
            # Inject restrt inside $name
            printf "        restrt = \047RESTRT\047\n"
            next 
        }

        # Reset block flags
        /\$end/ { in_para=0; in_name=0; print; next }

        # 5. Handle indexed file paths (skip comments)
        /boundn|saltn|tempn|vel3dn|wind|rain|qflux|saltin|tempin/ {
            if ($0 ~ /\.(dat|fem|txt)/ && $0 !~ /^[ \t]*(!|#)/) {
                if (match($0, /([\047][^\047]+[\047]|"[^"]+")/)) {
                    full_val_quoted = substr($0, RSTART, RLENGTH)
                    full_val = substr($0, RSTART+1, RLENGTH-2)
                    n = split(full_val, path_parts, "/")
                    fname = path_parts[n]
                    
                    base_n = fname; sub(/\.(dat|fem|txt)$/, "", base_n)
                    ext = (fname ~ /\.dat$/) ? ".dat" : (fname ~ /\.fem$/ ? ".fem" : ".txt")
                    indexed = base_n "_" idx ext
                    
                    if (system("[ -f " da_dir "/" indexed " ]") == 0) {
                        sub(full_val_quoted, "\047" indexed "\047")
                    } else if (system("[ -f " da_dir "/" fname " ]") == 0) {
                        sub(full_val_quoted, "\047" fname "\047")
                    }
                }
            }
        }

        # Print all other lines
        { print }
        ' "$base_str" > "$DA_DIR/$skel_file"
        
        [ -f "$rst_file" ] && mv "$rst_file" "$DA_DIR/"
        echo "$skel_file $rst_file" >> "$DA_DIR/ens_list.txt"
    done
}

#===================
handle_observations() {
    echo ">>> Step 5: Observations & Analysis Timeline"
    
    if [ -f "$OBS_DIR/obs_list.txt" ]; then
        echo "Found obs_list.txt. Linking observation data files..."
        cp "$OBS_DIR/obs_list.txt" "$DA_DIR/"
    else
        echo "WARNING: obs_list.txt not found in $OBS_DIR."
        read -p "--- [INPUT] Would you like to create it now? (yes/no): " create_obs
        
        if [[ "$create_obs" == "yes" ]]; then
            echo "--- Interactive Creation of obs_list.txt ---"
            #echo "# OBS_TYPE FILENAME X Y Z STD RHO" > "$DA_DIR/obs_list.txt"
            
            while true; do
                echo "Enter observation details (or type 'done' to finish):"
                read -p "Type (e.g., 0DSAL): " obs_type
                [[ "$obs_type" == "done" ]] && break
                
                read -p "Filename: " obs_file
                read -p "X (Lon): " obs_x
                read -p "Y (Lat): " obs_y
                read -p "Z (Depth): " obs_z
                read -p "STD (Std Dev): " obs_std
                read -p "RHO (Correlation Local Analysis): " obs_rho
                
                echo "$obs_type $obs_file $obs_x $obs_y $obs_z $obs_std $obs_rho" >> "$DA_DIR/obs_list.txt"
                echo "Row added."
            done
            echo "obs_list.txt created successfully in $DA_DIR."
        else
            echo "Data assimilation will not be possible without obs_list.txt."
        fi
    fi

    for f in "$OBS_DIR"/*.{dat,txt}; do
        [ -e "$f" ] && ln -sf "$f" "$DA_DIR/"
    done

    if [ -f "$ENKF_DIR/make_antime_list" ] && [ -f "$DA_DIR/obs_list.txt" ]; then
        echo "Configuring the Analysis list (frequency of DA steps)..."
        echo "Initial date: $R_DATE"
        read -p "--- [INPUT] Enter simulation final date (yyyy-mm-dd::HH:MM:SS): " T_END
        read -p "--- [INPUT] Enter minimum delta-time between analyses [seconds]: " T_DT
        (cd "$DA_DIR" && "$ENKF_DIR/make_antime_list" "$R_DATE" "$T_END" "$T_DT" "obs_list.txt")
    fi
}

#===================
run_enkf() {
    echo ">>> Step 6: Final Execution Settings"
    echo "Method 11-13: Stochastic EnKF | 21-23: Deterministic/Square-Root EnKF"
    read -p "--- [INPUT] Choose DA Method (e.g., 21): " DA_METHOD
    read -p "--- [INPUT] Apply Localization? (0=No, 1=Yes): " DA_LOC
    read -p "--- [INPUT] Number of parallel threads (nomp): " DA_THREADS
    read -p "--- [INPUT] Save all members output (required for enKS Smoothing)? (0/1): " DA_OUT
    read -p "--- [INPUT] Forecast execution engine (mpi|omp): " PMODE

    [ -f "$DET_DIR/lbound.dat" ] && cp "$DET_DIR/lbound.dat" "$DA_DIR/"
    
    read -p "--- [INPUT] Run the DA code? (y/n): " do_p
    if [ "$do_p" == "y" ]; then
       echo "Setup complete. Moving to $DA_DIR to launch the simulation."
       cd "$DA_DIR" 
       bash "$ENKF_DIR/enKF.sh" "$DA_METHOD" "$DA_LOC" "$DA_THREADS" "$DA_OUT" "$PMODE"
    else
       echo "Go to the DA dir and run: $ENKF_DIR/enKF.sh $DA_METHOD $DA_LOC $DA_THREADS $DA_OUT" "$PMODE"
    fi
}

# Main Execution Flow
check_env; ask_params; perturb_forcings_interactive; create_ensemble_data; handle_observations; run_enkf

