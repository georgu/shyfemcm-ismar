#!/usr/bin/env bash

# --- Path Deduction ---
ENKF_DIR=$(dirname "$(readlink -f "$0")")

# --- Usage Function ---
usage() {
    echo "Usage: $0 <fileheat> <nrens> <std_sol> <std_temp> <std_humi> <std_cloud> <Tau>"
    echo ""
    echo "Arguments:"
    echo "  fileheat  : Input file (Time Solar Rad Temp Humid Cloud)"
    echo "  nrens     : Number of ensemble members"
    echo "  std_sol   : STD for Solar Radiation"
    echo "  std_temp  : STD for Air Temperature"
    echo "  std_humi  : STD for Relative Humidity"
    echo "  std_cloud : STD for Cloud Cover"
    echo "  Tau       : Correlation time (Tau)"
    echo ""
    exit 1
}

# --- Check Arguments ---
if [ "$#" -ne 7 ]; then
    usage
fi

fileheat=$1
nrens=$2
s_sol=$3
s_temp=$4
s_humi=$5
s_cloud=$6
tau=$7

if [ ! -f "$fileheat" ]; then
    echo "ERROR: Input file '$fileheat' not found."
    exit 1
fi

fname="${fileheat##*/}"
bname="${fname%.*}"
ext="${fname##*.}"

# --- 1. Extract and Perturb each variable separately ---
# Column 2: Solar Radiation (Min: 0, Max: 1200)
awk '{print $1, $2}' "$fileheat" > sol.dat
"$ENKF_DIR/perturbe_ts" sol.dat "$nrens" "$s_sol" "$tau" 0. 1200.

# Column 3: Air Temperature (Min: -60, Max: 50)
awk '{print $1, $3}' "$fileheat" > temp.dat
"$ENKF_DIR/perturbe_ts" temp.dat "$nrens" "$s_temp" "$tau" -60. 50.

# Column 4: Relative Humidity (Min: 0, Max: 100)
awk '{print $1, $4}' "$fileheat" > humi.dat
"$ENKF_DIR/perturbe_ts" humi.dat "$nrens" "$s_humi" "$tau" 0. 100.

# Column 5: Cloud Cover (Min: 0, Max: 1)
awk '{print $1, $5}' "$fileheat" > cloud.dat
"$ENKF_DIR/perturbe_ts" cloud.dat "$nrens" "$s_cloud" "$tau" 0. 1.

# --- 2. Reconstruct Ensemble Members ---
echo "Reconstructing $nrens heat flux ensemble members..."
for ((i=0; i<nrens; i++)); do
    idx=$(printf "%03d" $i)
    
    # Merge all perturbed files for current index
    # paste aligns sol, temp, humi, and cloud files
    # awk picks the time from the first and the values from the others
    paste "sol_$idx.dat" "temp_$idx.dat" "humi_$idx.dat" "cloud_$idx.dat" | \
    awk '{print $1, $2, $4, $6, $8}' > "${bname}_$idx.$ext"
done

# --- 3. Cleanup ---
rm -f sol.dat temp.dat humi.dat cloud.dat
rm -f sol_[0-9][0-9][0-9].dat temp_[0-9][0-9][0-9].dat 
rm -f humi_[0-9][0-9][0-9].dat cloud_[0-9][0-9][0-9].dat 2>/dev/null

echo "Done. Ensemble files created: ${bname}_XXX.${ext}"
