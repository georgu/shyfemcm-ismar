#!/usr/bin/env bash

# --- Path Deduction ---
ENKF_DIR=$(dirname "$(readlink -f "$0")")

# --- Usage Function ---
usage() {
    echo "Usage: $0 <filewp> <nrens> <isws> <ispress> <STD_w> <STD_p> <Tau>"
    echo ""
    echo "Arguments:"
    echo "  filewp  : Path to the wind/pressure file"
    echo "  nrens   : Number of ensemble members to create"
    echo "  isws    : Wind format (0: u,v components | 1: speed,dir)"
    echo "  ispress : Pressure presence (0: no pressure | 1: with pressure)"
    echo "  STD_w   : Standard Deviation for wind speed perturbation"
    echo "  STD_p   : Standard Deviation for pressure perturbation"
    echo "  Tau     : Correlation time (Tau)"
    echo ""
    echo "Example: $0 wind.dat 10 0 1 1.5 100.0 86400"
    exit 1
}

# --- Check Arguments ---
if [ "$#" -ne 7 ]; then
    echo "ERROR: Wrong number of arguments."
    usage
fi

filewp=$1
nrens=$2
isws=$3
ispress=$4
STD_w=$5
STD_p=$6
Tau=$7

# Check if input file exists
if [ ! -f "$filewp" ]; then
    echo "ERROR: Input file '$filewp' not found."
    exit 1
fi

bname=$(echo "$filewp" | cut -d "." -f 1)
ext=$(echo "$filewp" | cut -d "." -f 2)

# --- 1. Extract Wind Speed (Magnitude) for perturbation ---
if [ "$isws" -eq "0" ]; then
    # Input is u, v components: calculate ws = sqrt(u^2 + v^2)
    awk '{print $1, sqrt($2^2 + $3^2)}' "$filewp" > ws.txt
else    
    # Input is already wind speed and direction
    awk '{print $1, $2}' "$filewp" > ws.txt
fi      

# Perturb the magnitude using the scalar perturbation tool
# Limits set to 0.0 and 60.0 m/s
"$ENKF_DIR/perturbe_ts" ws.txt "$nrens" "$STD_w" "$Tau" 0. 60.

# --- 2. Handle Pressure if required ---
if [ "$ispress" -eq "1" ]; then
   # Extract column 1 (time) and 4 (pressure)
   awk '{print $1, $4}' "$filewp" > press.txt
   # Perturb pressure (limits roughly 90000 to 110000 Pa)
   "$ENKF_DIR/perturbe_ts" press.txt "$nrens" "$STD_p" "$Tau" 90000. 110000.
fi 

# --- 3. Reconstruct Ensemble Members ---
echo "Reconstructing $nrens ensemble members..."
for ffile in ws_[0-9][0-9][0-9].txt; do
    idx=${ffile:3:3}
    
    if [ "$isws" -eq "0" ]; then
        # CASE U, V: Use the ratio (perturbed_ws / original_ws) as a multiplier
        # to scale components while preserving the original direction.
        paste "$filewp" "$ffile" | awk -v ispr="$ispress" '{
            orig_ws = sqrt($2^2 + $3^2);
            new_ws  = $5;
            ratio = (orig_ws > 0.001) ? new_ws / orig_ws : 0;
            
            if (ispr == "1") {
                # Intermediate step for pressure merge
                printf "%s %f %f ", $1, $2*ratio, $3*ratio;
            } else {
                # Final format (time, u, v)
                printf "%s %f %f\n", $1, $2*ratio, $3*ratio;
            }
        }' > tmp_u_v_$idx.txt

        if [ "$ispress" -eq "1" ]; then
            # Merge scaled u,v with perturbed pressure
            paste tmp_u_v_$idx.txt "press_$idx.txt" | awk '{print $1, $2, $3, $5}' > "${bname}_$idx.$ext"
        else
            mv tmp_u_v_$idx.txt "${bname}_$idx.$ext"
        fi
        rm -f tmp_u_v_$idx.txt

    else
        # CASE WS, DIR: Simply replace wind speed and keep original direction
        if [ "$ispress" -eq "1" ]; then
            paste "$filewp" "$ffile" "press_$idx.txt" | awk '{print $1, $5, $3, $7}' > "${bname}_$idx.$ext"
        else
            paste "$filewp" "$ffile" | awk '{print $1, $5, $3}' > "${bname}_$idx.$ext"
        fi
    fi
done

# --- 4. Cleanup temporary files ---
rm -f ws.txt press.txt ws_[0-9][0-9][0-9].txt press_[0-9][0-9][0-9].txt 2>/dev/null
echo "Done. Ensemble files created: ${bname}_XXX.${ext}"

