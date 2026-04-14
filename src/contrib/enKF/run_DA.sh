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
create_ensemble_data() {
    echo ">>> Step 3: Skeleton and Restart File Generation"
    local base_str=$(ls *.str 2>/dev/null | head -n 1)
    local base_rst=$(ls *.rst 2>/dev/null | head -n 1)
    local rst_name="${base_rst%.*}"

    local basin_name=$(awk '/\$title/{getline; getline; getline; print $1; exit}' "$base_str")
    local basin_f="${basin_name}.bas"
    [ -f "$DET_DIR/$basin_f" ] && cp -p "$DET_DIR/$basin_f" "$DA_DIR/"

    if [ -f "$base_rst" ]; then
        echo "Found base restart file: $base_rst"
        read -p "--- [INPUT] Do you want to PERTURB the initial state ($base_rst) to create ensemble spread? (y/n): " p_rst
        if [[ "$p_rst" == "y" ]]; then
            echo "Displaying info for $base_rst via rstinf..."
            "$SHYFEM_DIR/bin/rstinf" "$base_rst" 
            read -p "    > Analysis start date (FORMAT: yyyy-mm-dd::HH:MM:SS): " R_DATE
            read -p "    > Number of vertical levels in the model (nlv): " R_NLV
            read -p "    > Amplitude of SSH perturbation [meters] (errz): " R_ERRZ
            read -p "    > Amplitude of Temperature perturbation [°C] (errt): " R_ERRT
            read -p "    > Amplitude of Salinity perturbation [psu] (errs): " R_ERRS
            "$ENKF_DIR/perturbe_state" "$basin_f" "$base_rst" "$R_DATE" "$R_NLV" "$NRENS" "$R_ERRZ" "$R_ERRT" "$R_ERRS"
        else
            read -p "--- [INPUT] Enter the simulation start date (FORMAT: yyyy-mm-dd::HH:MM:SS): " R_DATE
            echo "Creating identical ensemble members from the deterministic restart..."
            for ((i=0; i<NRENS; i++)); do cp "$base_rst" "${rst_name}_$(printf "%03d" $i).rst"; done
        fi
    fi

    echo "Generating individual member Skeletons (.skel) from $base_str..."
    for ((i=0; i<NRENS; i++)); do
        idx=$(printf "%03d" $i)
        local skel_file="member_${idx}.skel"
        local rst_file="${rst_name}_${idx}.rst"

        awk -v idx="$idx" -v da_dir="$DA_DIR" '
        /\$title/ { print; t_step=1; next }
        t_step==1 { print "     Sim with DA"; t_step++; next }
        t_step==2 { print "     NAMESIM"; t_step++; next }
        t_step==3 { print $0; t_step=0; next }
        /\$para/ { print; print "        nocheck = 1"; next }
        /it(rst|mext|mout|mcon) *=/ { sub(/=.*/, "= \047ITANF\047"); print; next }
        /itanf  *=/ { print "        itanf  = \047ITANF\047   itend  = \047ITEND\047"; next }
        /idtrst *=/ { print "        idtrst = -1"; next }
        /restrt *=/ { print "        restrt = \047RESTRT\047"; next }
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
                        sub(full_val_quoted, "\047" indexed "\047"); print
                    } else if (system("[ -f " da_dir "/" fname " ]") == 0) {
                        sub(full_val_quoted, "\047" fname "\047"); print
                    } else { next }
                    next
                }
            }
        }
        { print }' "$base_str" > "$DA_DIR/$skel_file"
        
        [ -f "$rst_file" ] && mv "$rst_file" "$DA_DIR/"
        echo "$skel_file $rst_file" >> "$DA_DIR/ens_list.txt"
    done
}

#===================
perturb_forcings_interactive() {
    echo ">>> Step 4: Boundary & Forcing Perturbation"
    local base_str=$(ls *.str 2>/dev/null | head -n 1)
    local files=$(awk '/\$end/ {in_sec=0} /\$[a-zA-Z0-9]+/ {in_sec=1} in_sec && /boundn|saltn|tempn|vel3dn|wind|rain|qflux|saltin|tempin/ {if (match($0, /[\047"][^\047"]+\.(dat|fem|txt)[\047"]/)) {print substr($0, RSTART+1, RLENGTH-2)}}' "$base_str" | sort -u)

    for f in $files; do
        [ ! -f "$DET_DIR/$f" ] && continue
        local full_name=$(basename "$f")
        local name="${full_name%.*}"
        local ext="${f##*.}"

        echo "------------------------------------------------"
        echo "Detected forcing file: $full_name"
        read -p "--- [INPUT] Perturb this forcing to create ensemble variations? (y/n): " do_p
        if [[ "$do_p" == "y" ]]; then
            if [[ "$ext" == "dat" ]]; then
                read -p "      > Time Series Perturbation (STD MIN MAX Tau_rn): " TS_STD TS_MIN TS_MAX TS_TAU
                "$ENKF_DIR/perturbe_ts" "$f" "$NRENS" "$TS_STD" "$TS_TAU" "$TS_MIN" "$TS_MAX"
            elif [[ "$full_name" == *"wind"* || "$full_name" == *"press"* ]]; then
                read -p "      > Wind/Pressure Params (Method[1:P+W, 2:W] STD Tau_rn): " W_TYPE W_STD W_TAU
                "$ENKF_DIR/perturbe_fem_wind_press" "$f" "$NRENS" "$W_TYPE" "$W_STD" "$W_TAU"
            elif [[ "$full_name" == *"heat"* || "$full_name" == *"qflux"* ]]; then
                read -p "      > Heat Flux Params (STD_Solar STD_Temp STD_Humi STD_Cloud Tau_rn): " H_S H_A H_H H_C H_T
                "$ENKF_DIR/perturbe_fem_heat" "$f" "$NRENS" "$H_S" "$H_A" "$H_H" "$H_C" "$H_T"
            else
                read -p "      > Generic Scalar Params (STD MIN MAX Tau_rn): " S_STD S_VMIN S_VMAX S_TAU
                "$ENKF_DIR/perturbe_fem_scalar" "$f" "$NRENS" "$S_STD" "$S_VMIN" "$S_VMAX" "$S_TAU"
            fi
            mv ${name}_[0-9][0-9][0-9].${ext} "$DA_DIR/" 2>/dev/null
        else
            echo "      - Copying original $full_name as static forcing for all members."
            cp "$DET_DIR/$f" "$DA_DIR/$full_name"
        fi
    done
}

#===================
handle_observations() {
    echo ">>> Step 5: Observations & Analysis Timeline"
    if [ -f "$OBS_DIR/obs_list.txt" ]; then
        echo "Found obs_list.txt. Linking observation data files..."
        cp "$OBS_DIR/obs_list.txt" "$DA_DIR/"
        ln -sf "$OBS_DIR"/*.{dat,txt} "$DA_DIR/" 2>/dev/null
    else
        echo "WARNING: obs_list.txt not found in $OBS_DIR. Data assimilation will not be possible."
    fi

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

    [ -f "$DET_DIR/lbound.dat" ] && cp "$DET_DIR/lbound.dat" "$DA_DIR/"
    
    echo "Setup complete. Moving to $DA_DIR to launch the simulation."
    cd "$DA_DIR" && bash "$ENKF_DIR/enKF.sh" "$DA_METHOD" "$DA_LOC" "$DA_THREADS" "$DA_OUT"
}

# Main Execution Flow
check_env; ask_params; create_ensemble_data; perturb_forcings_interactive; handle_observations; run_enkf

