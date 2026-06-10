#!/usr/bin/env bash

# --- Path Deduction ---
# Locate the directory of the script and the perturbe_ts executable
ENKF_DIR=$(dirname "$(readlink -f "$0")")

# --- Usage Function ---
usage() {
    echo "Usage: $0 <filewp> <nrens> <isws> <ispress> <STD_w> <STD_p> <Tau>"
    echo ""
    echo "Arguments:"
    echo "  filewp  : Path to the wind/pressure file (Time format: YYYY-mm-dd::HH:MM:SS)"
    echo "  nrens   : Number of ensemble members to create"
    echo "  isws    : Wind format (0: u,v components | 1: speed,dir)"
    echo "  ispress : Pressure presence (0: no pressure | 1: with pressure)"
    echo "  STD_w   : Standard Deviation for wind speed perturbation"
    echo "  STD_p   : Standard Deviation for pressure perturbation"
    echo "  Tau     : Correlation time (Tau in seconds)"
    exit 1
}

if [ "$#" -ne 7 ]; then usage; fi

filewp=$1
nrens=$2
isws=$3
ispress=$4
STD_w=$5
STD_p=$6
Tau=$7

# Check if executable exists
if [ ! -x "$ENKF_DIR/perturbe_ts" ]; then
    echo "ERROR: Executable 'perturbe_ts' not found in $ENKF_DIR"
    exit 1
fi

fname="${filewp##*/}"
bname="${fname%.*}"
ext="${fname##*.}"

# --- 1. Extract Wind Speed (Magnitude) ---
if [ "$isws" -eq "0" ]; then
    # Input is Time U V [P]
    awk '{print $1, sqrt($2^2 + $3^2)}' "$filewp" > ws.dat
else    
    # Input is Time WS DIR [P]
    awk '{print $1, $2}' "$filewp" > ws.dat
fi      

# Perturb the magnitude (vmin=0.0 as requested, vmax=60.0)
"$ENKF_DIR/perturbe_ts" ws.dat "$nrens" "$STD_w" "$Tau" 0.0 60.0

# --- 2. Handle Pressure if required ---
if [ "$ispress" -eq "1" ]; then
   awk '{print $1, $4}' "$filewp" > press.dat
   "$ENKF_DIR/perturbe_ts" press.dat "$nrens" "$STD_p" "$Tau" 90000. 110000.
fi 

# --- 3. Reconstruct Ensemble Members ---
echo "Reconstructing $nrens ensemble members..."
for ffile in ws_[0-9][0-9][0-9].dat; do
    idx=${ffile:3:3}
    
    if [ "$isws" -eq "0" ]; then
        # CASE U, V: Ratio Scaling with Floor Clipping
        # Paste: $1=Time, $2=U_orig, $3=V_orig, $4=P_orig, $5=Time_rep, $6=WS_pert
        paste "$filewp" "$ffile" | awk -v ispr="$ispress" '{
            u_orig = $2; v_orig = $3;
            orig_ws = sqrt(u_orig^2 + v_orig^2);
            new_ws  = $6; 
            
            # Applying a floor of 0.5 m/s to avoid null components when STD > signal
            if (new_ws < 0.5) new_ws = 0.5;

            # Safe scaling (orig_ws > 0 as it comes from observations)
            ratio = new_ws / orig_ws;
            u_new = u_orig * ratio;
            v_new = v_orig * ratio;

            if (ispr == "1") printf "%s %12.5f %12.5f %12.5f\n", $1, u_new, v_new, $4;
            else printf "%s %12.5f %12.5f\n", $1, u_new, v_new;
        }' > tmp_scaled_$idx.dat

        if [ "$ispress" -eq "1" ]; then
            # Replace original pressure with perturbed pressure from press_XXX.dat ($6)
            paste tmp_scaled_$idx.dat "press_$idx.dat" | awk '{print $1, $2, $3, $6}' > "${bname}_$idx.$ext"
        else
            mv tmp_scaled_$idx.dat "${bname}_$idx.$ext"
        fi
        rm -f tmp_scaled_$idx.dat
    else
        # CASE WS, DIR: Replace WS (with 0.5 floor), keep original DIR
        # Paste: $1=Time, $2=WS_orig, $3=DIR_orig, $4=P_orig, $5=Time_rep, $6=WS_pert
        if [ "$ispress" -eq "1" ]; then
            paste "$filewp" "$ffile" "press_$idx.dat" | awk '{
                new_ws = ($6 < 0.5) ? 0.5 : $6;
                print $1, new_ws, $3, $8
            }' > "${bname}_$idx.$ext"
        else
            paste "$filewp" "$ffile" | awk '{
                new_ws = ($6 < 0.5) ? 0.5 : $6;
                print $1, new_ws, $3
            }' > "${bname}_$idx.$ext"
        fi
    fi
done

# --- 4. Cleanup ---
rm -f ws.dat press.dat ws_[0-9][0-9][0-9].dat press_[0-9][0-9][0-9].dat 2>/dev/null
echo "Ensemble generation complete: ${bname}_XXX.${ext}"

