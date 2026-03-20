#!/usr/bin/env bash
#------------------------------------------------------------------------------
# PARALLEL SHYFEM ENSEMBLE PROCESSOR
# 1. SPLIT: Parallel execution of shyelab in dedicated temp folders.
# 2. CLEAN: Remove headers (#) and rename fragments to 1-based stations.
# 3. MERGE: Temporal join of an* files for each member/station.
# 4. STATS: Cross-ensemble statistics (Mean, Std, Max, Min) per station.
#------------------------------------------------------------------------------

# --- CONFIGURATION ---
MAX_JOBS=10           # Number of parallel ensemble members to process
SCRIPTPATH="$(cd "$(dirname "$0")" && pwd)"
FEMDIR="${SCRIPTPATH}/../../.."
SHYELAB="${FEMDIR}/bin/shyelab"

# Check dependencies
[ -x "$SHYELAB" ] || { echo "FATAL: shyelab not found at $SHYELAB"; exit 1; }

# Variables to look for
vars='zeta.2d velx.2d vely.2d velx.3d vely.3d speed.2d speed.3d dir.2d dir.3d all.2d temp.2d temp.3d salt.2d salt.3d'

# 1. Identify unique ensemble members from filenames (e.g., en00000b)
members=$(ls an*_en*.ext 2>/dev/null | sed 's/.*_\(en.*\)\.ext/\1/' | sort -u)

# Function to process a single ensemble member
process_member() {
    local mem=$1
    local tmp_dir="tmp_${mem}"
    mkdir -p "$tmp_dir"

    echo "[Member $mem] Starting parallel split..."

    # Process each temporal file (.ext) for this member
    for efile in an*_"${mem}".ext; do
        [ -f "$efile" ] || continue
        basefile=$(basename "$efile" .ext)

        # Create a symbolic link inside the temp folder to avoid moving large files
        ln -sf "../$efile" "$tmp_dir/$efile"

        # Run shyelab split inside the isolated temp directory
        (cd "$tmp_dir" && "$SHYELAB" -split "$efile" > /dev/null 2>&1)

        # Clean headers and rename fragments to 1-based index
        for vv in $vars; do
            if ls "$tmp_dir"/${vv}.[0-9]* >/dev/null 2>&1; then
                for fl in "$tmp_dir"/${vv}.[0-9]*; do
                    # Extract 0-based index and convert to 1-based station
                    idx=$(echo "$fl" | rev | cut -d'.' -f1 | rev)
                    new_idx=$((idx + 1))

                    # Remove # headers and save as temporal station file
                    sed '/^#/d' "$fl" > "${basefile}_${vv}_st${new_idx}.ts"
                    rm -f "$fl"
                done
            fi
        done
        rm -f "$tmp_dir/$efile" # Remove symlink
    done
    rmdir "$tmp_dir" # Cleanup temp folder

    # Phase 2: Immediate temporal merge for this member
    echo "[Member $mem] Merging temporal fragments..."
    active_vars=$(ls an*_"${mem}"_*_st*.ts 2>/dev/null | sed "s/an[0-9]*_${mem}_\(.*\)_st[0-9]*.ts/\1/" | sort -u)

    for vv in $active_vars; do
        stations=$(ls an*_"${mem}"_"${vv}"_st*.ts 2>/dev/null | sed "s/.*_st\([0-9]*\)\.ts/\1/" | sort -un)
        for st in $stations; do
            output="TOTAL_${mem}_${vv}_st${st}.ts"
            # Concatenate fragments in natural order (an001 -> an010)
            ls an*_"${mem}"_"${vv}"_st${st}.ts 2>/dev/null | sort -V | xargs cat > "$output"
            # Cleanup fragments after successful merge
            [ -s "$output" ] && rm -f an*_"${mem}"_"${vv}"_st${st}.ts
        done
    done
}

# --- PARALLEL EXECUTION ENGINE ---
count=0
for mem in $members; do
    process_member "$mem" &
    ((count++))

    # Limit number of background jobs
    if [ "$count" -ge "$MAX_JOBS" ]; then
        wait
        count=0
    fi
done
wait # Ensure all members are finished before calculating global stats

# --- PHASE 3: GLOBAL ENSEMBLE STATISTICS (Row-wise) ---
echo "--- Calculating Global Ensemble Statistics ---"
# Find all unique variable+station combinations
all_targets=$(ls TOTAL_*_st*.ts 2>/dev/null | sed 's/TOTAL_en[^_]*_\(.*\)\.ts/\1/' | sort -u)

for target in $all_targets; do
    # Collect all member files for this station/variable
    member_files=$(ls TOTAL_*_"${target}".ts 2>/dev/null | sort -V)
    n_files=$(echo "$member_files" | wc -w)

    if [ "$n_files" -gt 0 ]; then
        output_stats="STATS_ENSEMBLE_${target}.ts"
        echo "Generating $output_stats ($n_files members)"

        # Paste all members horizontally and compute stats per row
        paste $member_files | awk -v n="$n_files" '
        {
            sum=0; sumsq=0; max=-1e30; min=1e30;
            timestamp=$1;
            for (i=0; i<n; i++) {
                val=$( (i*2)+2 ); # Value is in the 2nd column of each member block
                sum+=val; sumsq+=val*val;
                if(val>max) max=val; if(val<min) min=val;
            }
            mean=sum/n;
            var=(sumsq/n)-(mean*mean);
            std=(var>0)?sqrt(var):0;
            # Output format: Time Mean StdDev Max Min
            printf "%s\t%.6f\t%.6f\t%.6f\t%.6f\n", timestamp, mean, std, max, min;
        }' > "$output_stats"
    fi
done

echo "Workflow completed successfully."
