#!/bin/bash

# --- CHECK INPUT ARGUMENT ---
if [ -z "$1" ]; then
    echo "Usage: $0 <base_filename>"
    echo "Example: $0 pe_parameters"
    exit 1
fi

BASE_NAME="$1"
OUT_DIR="./pe_time_series"
mkdir -p "$OUT_DIR"

echo "Checking available files matching base pattern: ${BASE_NAME}_anXXXXX_enXXXXX..."

# --- PHASE 1: AUTOMATICALLY DETECT AVAILABLE MEMBERS AND TIMESTEPS ---
# Extract unique 5-digit analysis steps (anXXXXX) based on the input prefix
AN_STEPS=$(ls "${BASE_NAME}"_an[0-9][0-9][0-9][0-9][0-9]_en[0-9][0-9][0-9][0-9][0-9]*.dat 2>/dev/null | \
           grep -o '_an[0-9]\{5\}' | sort -u | sed 's/_an//')

# Extract unique 5-digit ensemble members (enXXXXX) based on the input prefix
EN_MEMBERS=$(ls "${BASE_NAME}"_an[0-9][0-9][0-9][0-9][0-9]_en[0-9][0-9][0-9][0-9][0-9]*.dat 2>/dev/null | \
             grep -o '_en[0-9]\{5\}' | sort -u | sed 's/_en//')

if [ -z "$AN_STEPS" ] || [ -z "$EN_MEMBERS" ]; then
    echo "Error: No matching files found for pattern: ${BASE_NAME}_anXXXXX_enXXXXX*.dat"
    exit 1
fi

AN_ARRAY=($AN_STEPS)
EN_ARRAY=($EN_MEMBERS)
NUM_ENS=${#EN_ARRAY[@]}

echo "Found ${#AN_ARRAY[@]} timesteps and $NUM_ENS ensemble members."

# --- PHASE 2: EXTRACT PARAMETER LABELS & SUFFIX ---
FIRST_FILE=$(ls "${BASE_NAME}"_an[0-9][0-9][0-9][0-9][0-9]_en[0-9][0-9][0-9][0-9][0-9]*.dat 2>/dev/null | head -n 1)
SUFFIX=$(echo "$FIRST_FILE" | grep -o 'en[0-9]\{5\}.*' | sed 's/en[0-9]\{5\}//')

PARAMS=$(awk '{print $1}' "$FIRST_FILE")

echo "Suffix detected: $SUFFIX"
echo "Initializing output files for each parameter..."

# --- PHASE 3: INITIALIZE HEADERS FOR EACH PARAMETER FILE ---
for param in $PARAMS; do
    printf "# STEP\tMEAN\tSTD" > "$OUT_DIR/${param}_ts.dat"
    for en in "${EN_ARRAY[@]}"; do
        printf "\ten%s" "$en" >> "$OUT_DIR/${param}_ts.dat"
    done
    printf "\n" >> "$OUT_DIR/${param}_ts.dat"
done

# --- PHASE 4: PROCESSING TIME STEPS & COMPUTING STATS ---
echo "Processing time series..."

for an in "${AN_ARRAY[@]}"; do
    echo "  -> Processing analysis step: $an"
    
    # We use a standard awk logic tracking FILENAME switches to increment file index (file_idx)
    awk -v an_step="$an" -v nens="$NUM_ENS" -v out_dir="$OUT_DIR" '
    BEGIN {
        row = 0
        file_idx = 0
        last_filename = ""
    }
    {
        # Increment file index when awk switches to the next input file
        if (FILENAME != last_filename) {
            file_idx++
            last_filename = FILENAME
            row = 0
        }
        
        row++
        
        # Store label from the very first file processed
        if (file_idx == 1) {
           labels[row] = $1
        }
        
        # Save value into matrix
        values[row, file_idx] = $2
    }
    END {
        total_rows = row
        for (r = 1; r <= total_rows; r++) {
            sum = 0
            sq_sum = 0
            
            # Loop through the files collected
            for (m = 1; m <= nens; m++) {
                v = values[r, m]
                sum += v
                sq_sum += v * v
            }
            mean = sum / nens
            
            variance = (sq_sum / nens) - (mean * mean)
            if (variance < 0) variance = 0
            std = sqrt(variance)
            
            out_file = out_dir "/" labels[r] "_ts.dat"
            printf "%s\t%.6f\t%.6f", an_step, mean, std >> out_file
            for (m = 1; m <= nens; m++) {
                printf "\t%.6f", values[r, m] >> out_file
            }
            printf "\n" >> out_file
        }
    }
    ' "${BASE_NAME}"_an"${an}"_en[0-9][0-9][0-9][0-9][0-9]"${SUFFIX}"
done

echo "Done! All time series files are available in: $OUT_DIR"

