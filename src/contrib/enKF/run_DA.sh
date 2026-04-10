#!/usr/bin/env bash

# --- Path Deduction ---
ENKF_DIR=$(dirname "$(readlink -f "$0")")
SHYFEM_DIR=$(cd "$ENKF_DIR/../../.." && pwd)
DET_DIR=$(pwd)
DA_DIR="$DET_DIR/sim_DA"
OBS_DIR="$DET_DIR/obs"

check_env() {
    echo ">>> Step 1: Environment Setup"
    rm -fr "$DA_DIR"
    mkdir -p "$DA_DIR"
    [ ! -f "$ENKF_DIR/enKF.sh" ] && { echo "Error: enKF.sh missing."; exit 1; }
}

ask_params() {
    echo ">>> Step 2: Main Parameters"
    read -p "Enter number of ensemble members (NRENS): " NRENS
    [ $((NRENS % 2)) -eq 0 ] && NRENS=$((NRENS + 1)) && echo "Adjusted to $NRENS (odd)."

    echo "DA Execution Settings:"
    read -p "Method (11-23): " DA_METHOD
    read -p "Localisation (0/1): " DA_LOC
    read -p "Threads: " DA_THREADS
    read -p "Output (0/1): " DA_OUT
}

create_ensemble_data() {
    echo ">>> Step 3: Skeleton and Restart Generation"
    local base_str=$(ls *.str 2>/dev/null | head -n 1)
    local base_rst=$(ls *.rst 2>/dev/null | head -n 1)
    local rst_name="${base_rst%.*}"

    local basin_name=$(awk '/\$title/{getline; getline; getline; print $1; exit}' "$base_str")
    local basin_f="${basin_name}.bas"

    if [ -f "$base_rst" ]; then
        read -p "Perturb initial state $base_rst (NRENS=$NRENS)? (y/n): " p_rst
        if [[ "$p_rst" == "y" ]]; then
            read -p "Date (yyyy-mm-dd::HH:MM:SS): " R_DATE
            read -p "Vertical levels (nlv): " R_NLV
            read -p "Amp SSH (errz): " R_ERRZ
            read -p "Amp Temp (errt): " R_ERRT
            read -p "Amp Salt (errs): " R_ERRS
            "$ENKF_DIR/perturbe_state" "$basin_f" "$base_rst" "$R_DATE" "$R_NLV" "$NRENS" "$R_ERRZ" "$R_ERRT" "$R_ERRS"
        else
            for ((i=0; i<NRENS; i++)); do cp "$base_rst" "${rst_name}_$(printf "%03d" $i).rst"; done
        fi
    fi

    > "$DA_DIR/ens_list.txt"
    for ((i=0; i<NRENS; i++)); do
        idx=$(printf "%03d" $i)
        local skel_file="member_${idx}.skel"
        local rst_file="${rst_name}_${idx}.rst"

        awk -v idx="$idx" '
        /\$title/ { print; print "     Sim with DA"; print "     NAMESIM"; getline; getline; next }
        /itrst =/ || /itmext =/ || /itmout =/ || /itmcon =/ { sub(/=.*/, "= \047ITANF\047"); print; next }
        /itanf  =/ { print "        itanf  = \047ITANF\047   itend  = \047ITEND\047"; next }
        /idtrst =/ { print "        idtrst = -1"; next }
        /it[lr][ae]nf =/ || /itmlgr =/ { gsub(/\047[^\047]+\047/, "\047ITANF\047"); print; next }
        /restrt =/ { print "        restrt = \047RESTRT\047"; next }
        /boundn|saltn|tempn|vel3dn|wind|rain|qflux/ {
            if ($0 ~ /\.(dat|fem)/ && $0 !~ /^[ \t]*(!|#)/) {
                sub(/\.dat/, "_" idx ".dat"); sub(/\.fem/, "_" idx ".fem")
            }
        }
        { print }' "$base_str" > "$DA_DIR/$skel_file"

        [ -f "$rst_file" ] && mv "$rst_file" "$DA_DIR/"
        echo "$skel_file $rst_file" >> "$DA_DIR/ens_list.txt"
    done
}

perturb_forcings_interactive() {
    echo ">>> Step 4: Interactive Forcing Perturbation (Batch Mode)"
    local base_str=$(ls *.str 2>/dev/null | head -n 1)
    local files=$(grep -E '(boundn|saltn|tempn|vel3dn|wind|rain|qflux)' "$base_str" | \
                  grep -vE '^[ \t]*(!|#|\$)' | sed "s/.*= *['\"]//;s/['\"].*//;s|.*/||" | sort -u)

    for f in $files; do
        [ ! -f "$DET_DIR/$f" ] && continue
        local ext="${f##*.}"
        local name="${f%.*}"

        read -p "Perturb forcing $f? (y/n): " do_p
        if [[ "$do_p" == "y" ]]; then
            if [[ "$ext" == "dat" ]]; then
                read -p "TS Params (STD Tau MIN MAX): " TS_STD TS_TAU TS_MIN TS_MAX
                "$ENKF_DIR/perturbe_ts" "$f" "$NRENS" "$TS_STD" "$TS_TAU" "$TS_MIN" "$TS_MAX"
            elif [[ "$f" == *"wind"*.fem || "$f" == *"press"*.fem  || "$f" == *"meteo"*.fem ]]; then
                read -p "Wind Params (Type[1:P+W, 2:W] STD Tau): " W_TYPE W_STD W_TAU
                "$ENKF_DIR/perturbe_fem_wind_press" "$f" "$NRENS" "$W_TYPE" "$W_STD" "$W_TAU"
            elif [[ "$f" == *"heat"*.fem || "$f" == *"qflux"*.fem ]]; then
                read -p "Heat Params (STD_solar STD_temp STD_humi STD_cloud Tau): " H_S H_A H_H H_C H_T
                "$ENKF_DIR/perturbe_fem_heat" "$f" "$NRENS" "$H_S" "$H_A" "$H_H" "$H_C" "$H_T"
            else
                read -p "Scalar Params (STD_r vmin vmax Tau): " S_STD S_VMIN S_VMAX S_TAU
                "$ENKF_DIR/perturbe_fem_scalar" "$f" "$NRENS" "$S_STD" "$S_VMIN" "$S_VMAX" "$S_TAU"
            fi
        else
            for ((i=0; i<NRENS; i++)); do cp "$f" "${name}_$(printf "%03d" $i).${ext}"; done
        fi
        mv ${name}_[0-9][0-9][0-9].${ext} "$DA_DIR/" 2>/dev/null
    done
}

run_enkf() {
    echo ">>> Step 5: Finalizing and Launch"
    [ -d "$OBS_DIR" ] && ln -sf "$OBS_DIR"/* "$DA_DIR/"
    [ -f "$DET_DIR/lbound.dat" ] && cp "$DET_DIR/lbound.dat" "$DA_DIR/"
    
    cd "$DA_DIR" || exit
    [ -f "$ENKF_DIR/make_antime_list" ] && "$ENKF_DIR/make_antime_list"
    
    read -p "Ready. Launch enKF.sh? (y/n): " start_sim
    [[ "$start_sim" == "y" ]] && bash "$ENKF_DIR/enKF.sh" "$DA_METHOD" "$DA_LOC" "$DA_THREADS" "$DA_OUT"
}

check_env; ask_params; create_ensemble_data; perturb_forcings_interactive; run_enkf

