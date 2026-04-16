#!/usr/bin/env bash

# --- Path Deduction ---
# Locate the directory where the script and the perturbe_ts executable reside
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

# Check if the executable exists
if [ ! -x "$ENKF_DIR/perturbe_ts" ]; then
    echo "ERROR: Executable 'perturbe_ts' not found or not executable in $ENKF_DIR"
    exit 1
fi

# Extract filename parts for output naming
fname="${filewp##*/}"
bname="${fname%.*}"
ext="${fname##*.}"

# --- 1. Extract Wind Speed (Magnitude) for perturbation ---
if [ "$isws" -eq "0" ]; then
    # Input is u, v: calculate ws = sqrt(u^2 + v^2). Time is $1.
    awk '{print $1, sqrt($2^2 + $3^2)}' "$filewp" > ws.dat
else
    # Input is already wind speed and direction. Extract Time and WS.
    awk '{print $1, $2}' "$filewp" > ws.dat
fi

# Perturb the magnitude using the scalar perturbation tool
# Limits for wind speed are set between 0.0 and 60.0 m/s
"$ENKF_DIR/perturbe_ts" ws.dat "$nrens" "$STD_w" "$Tau" 0. 60.

# --- 2. Handle Pressure if required ---
if [ "$ispress" -eq "1" ]; then
   # Extract column 1 (time) and 4 (pressure)
   awk '{print $1, $4}' "$filewp" > press.dat
   # Perturb pressure (physical limits roughly 90000 to 110000 Pa)
   "$ENKF_DIR/perturbe_ts" press.dat "$nrens" "$STD_p" "$Tau" 90000. 110000.
fi

# --- 3. Reconstruct Ensemble Members ---
echo "Reconstructing $nrens ensemble members..."
for ffile in ws_[0-9][0-9][0-9].dat; do
    # Extract the 3-digit index (e.g., 001 from ws_001.dat)
    idx=${ffile:3:3}

    if [ "$isws" -eq "0" ]; then
        # CASE U, V: Scale components while handling calm winds
        # After paste: $1=Time, $2=U_orig, $3=V_orig, $4=P_orig, $5=Time_rep, $6=WS_pert
        paste "$filewp" "$ffile" | awk -v ispr="$ispress" 'BEGIN{srand()} {
            u_orig = $2; v_orig = $3;
            orig_ws = sqrt(u_orig^2 + v_orig^2);
            new_ws  = $6;

            if (orig_ws > 0.1) {
                # Normal case: scale existing components (preserves direction)
                ratio = new_ws / orig_ws;
                u_new = u_orig * ratio;
                v_new = v_orig * ratio;
            } else if (new_ws > 0.1) {
                # Calm wind case: assign a random direction to the new wind magnitude
                angle = rand() * 2 * 3.14159265;
                u_new = new_ws * cos(angle);
                v_new = new_ws * sin(angle);
            } else {
                # Both original and perturbed are effectively zero
                u_new = 0.0; v_new = 0.0;
            }

            if (ispr == "1") {
                # Print temporary line with original pressure (updated in next step)
                printf "%s %f %f %f\n", $1, u_new, v_new, $4;
            } else {
                printf "%s %f %f\n", $1, u_new, v_new;
            }
        }' > tmp_scaled_$idx.dat

        if [ "$ispress" -eq "1" ]; then
            # Merge scaled u,v with perturbed pressure from press_XXX.dat
            # After paste: $1..$4 = data from tmp, $6 = perturbed pressure
            paste tmp_scaled_$idx.dat "press_$idx.dat" | awk '{print $1, $2, $3, $6}' > "${bname}_$idx.$ext"
        else
            mv tmp_scaled_$idx.dat "${bname}_$idx.$ext"
        fi
        rm -f tmp_scaled_$idx.dat

    else
        # CASE WS, DIR: Simply replace wind speed and keep original direction
        # After paste: $1=Time, $2=WS_orig, $3=DIR_orig, $4=P_orig, $5=Time, $6=WS_pert
        if [ "$ispress" -eq "1" ]; then
            paste "$filewp" "$ffile" "press_$idx.dat" | awk '{print $1, $6, $3, $8}' > "${bname}_$idx.$ext"
        else
            paste "$filewp" "$ffile" | awk '{print $1, $6, $3}' > "${bname}_$idx.$ext"
        fi
    fi
done

# --- 4. Cleanup temporary files ---
rm -f ws.dat press.dat ws_[0-9][0-9][0-9].dat press_[0-9][0-9][0-9].dat 2>/dev/null
echo "Done. Ensemble files created: ${bname}_XXX.${ext}"
