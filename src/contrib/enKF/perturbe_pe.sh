#!/bin/bash

# Ensure both the input file and NRENS arguments were provided on the command line
if [ -z "$1" ] || [ -z "$2" ]; then
    echo "Error: Missing arguments. Provide the input file and NRENS (e.g., ./script.sh file.dat 10)." >&2
    exit 1
fi

# Assign the command-line arguments to variables
INPUT_FILE="$1"
NRENS="$2"

# Verify that the specified source file actually exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file '$INPUT_FILE' not found!" >&2
    exit 1
fi

INPUT_FILE_B=$(basename $INPUT_FILE .dat)

echo "Starting ensemble generation for $NRENS members with 3-Sigma and Min/Max limits..."

# Initialize a dynamic seed for the random number generator using system time and Process ID
RANDOM_SEED=$(( ( $(date +%s) + $$ ) % 32768 ))

# --- MAIN LOOP ---
for ((i=0; i<NRENS; i++)); do
    # Format the member ID with leading zeros
    MEMBER_ID=$(printf "%05d" "$i")
    OUTPUT_FILE="${INPUT_FILE_B}_an00001_en${MEMBER_ID}.dat"
    
    echo "Generating: $OUTPUT_FILE"
    
    # Process the input file using AWK to handle 5 columns and mathematical filtering
    awk -v member="$i" -v seed="$RANDOM_SEED" '
    BEGIN {
        # Seed the random number generator uniquely for each ensemble member loop
        srand(seed + member); 
    }
    {
        # Skip empty lines or lines starting with a comment character (#)
        if ($1 == "" || $1 ~ /^#/) {
            print $0;
            next;
        }

        # Extract values from the 5 columns of the input file
        name = $1;
        value = $2;
        std_val = $3;
        min_val = $4;
        max_val = $5;

        # Check if this is the control member (member 000)
        if (member == 0) {
            # Member 000 remains completely unperturbed
            new_value = value;
        } else {
            # Box-Muller transform to generate a standard Gaussian random number (mean=0, std=1)
            u1 = rand();
            u2 = rand();
            
            # Prevent a division by zero / math error if u1 happens to be exactly 0
            if (u1 == 0) u1 = 0.00001; 
            
            # Generate the Gaussian noise using the cosine component of the Box-Muller formula
            g_noise = sqrt(-2.0 * log(u1)) * cos(2.0 * 3.141592653589793 * u2);
            
            # Calculate 3-sigma boundaries based on the current parameter
            sigma3_lower = value - (3.0 * std_val);
            sigma3_upper = value + (3.0 * std_val);

            # Apply the target formula: Value + (Standard_Gaussian_Noise * Column_3_STD)
            new_value = value + (g_noise * std_val);

            # --- 1. APPLY 3-SIGMA BOUNDARY LIMITS ---
            if (new_value < sigma3_lower) new_value = sigma3_lower;
            if (new_value > sigma3_upper) new_value = sigma3_upper;

            # --- 2. APPLY MIN/MAX HARD BOUNDARY LIMITS ---
            if (new_value < min_val) new_value = min_val;
            if (new_value > max_val) new_value = max_val;
        }

        # Output all 5 columns to maintain full consistency with the new input file structure
        printf "%s\t%.5f\t%s\t%s\t%s\n", name, new_value, std_val, min_val, max_val;
    }' "$INPUT_FILE" > "$OUTPUT_FILE"

    # --- AUTOMATED SKELSTR CALL ---
    # SkelStr "val1" "val2" "val3" "val4" "input_template.txt" "$OUTPUT_FILE" "final_output_${MEMBER_ID}.txt"

done

echo "Ensemble generation complete!"

