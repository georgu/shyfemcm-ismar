#!/usr/bin/env bash
#------------------------------------------------------------------------------
# PARALLEL SHYFEM ENSEMBLE PROCESSOR
#------------------------------------------------------------------------------

# --- CONFIGURATION ---
MAX_JOBS=10
SCRIPTPATH="$(cd "$(dirname "$0")" && pwd)"
FEMDIR="${SCRIPTPATH}/../../.."
SHYELAB="${FEMDIR}/bin/shyelab"

# --- HELP FUNCTION ---
show_help() {
    echo "Usage: $0 [MODE_PREFIX]"
    echo ""
    echo "Arguments:"
    echo "  MODE_PREFIX    The file prefix to process (e.g., 'an' or 'exp_name')."
    echo ""
    echo "Examples:"
    echo "  1. AN Mode (fragmented): $0 an"
    echo "     Merges multiple an*_en*b.ext fragments per member."
    echo ""
    echo "  2. Direct Mode:          $0 experiment"
    echo "     Each experiment_*.ext is already a full member."
    echo ""
    exit 0
}

if [ -z "$1" ]; then show_help; fi

# --- INPUT HANDLING ---
INPUT_NAME=$1
[ -x "$SHYELAB" ] || { echo "FATAL: shyelab not found at $SHYELAB"; exit 1; }

vars='zeta.2d velx.2d vely.2d velx.3d vely.3d speed.2d speed.3d dir.2d dir.3d all.2d temp.2d temp.3d salt.2d salt.3d'

if [ "$INPUT_NAME" == "an" ]; then
    echo "Mode: Temporal Fragments (an*_en*.ext)"
    MODE="AN"
    members=$(ls an*_en*.ext 2>/dev/null | sed 's/.*_\(en.*\)\.ext/\1/' | sort -u)
else
    echo "Mode: Direct Members (${INPUT_NAME}_*.ext)"
    MODE="NAME"
    members=$(ls ${INPUT_NAME}_*.ext 2>/dev/null | sed "s/${INPUT_NAME}_\(.*\)\.ext/\1/" | sort -u)
fi

[ -z "$members" ] && { echo "ERROR: No files found for prefix '$INPUT_NAME'."; exit 1; }

# --- CORE PROCESSING FUNCTION ---
process_member() {
    local mem=$1
    local tmp_dir="tmp_${mem}"
    mkdir -p "$tmp_dir"

    echo "[Member $mem] Splitting..."

    if [ "$MODE" == "AN" ]; then
        files_to_proc=$(ls an*_"${mem}".ext 2>/dev/null)
    else
        files_to_proc="${INPUT_NAME}_${mem}.ext"
    fi

    for efile in $files_to_proc; do
        [ -f "$efile" ] || continue
        ln -sf "../$efile" "$tmp_dir/$efile"

        # Run split
        (cd "$tmp_dir" && "$SHYELAB" -split "$efile" > /dev/null 2>&1)

        # Clean headers and move to processing area
        for vv in $vars; do
            if ls "$tmp_dir"/${vv}.[0-9]* >/dev/null 2>&1; then
                for fl in "$tmp_dir"/${vv}.[0-9]*; do
                    idx=$(echo "$fl" | rev | cut -d'.' -f1 | rev)
                    new_idx=$((idx + 1))
                    # Use a unique temporary name including the original file name to avoid collisions in AN mode
                    sed '/^#/d' "$fl" > "PROC_${mem}_${efile}_${vv}_st${new_idx}.ts"
                done
            fi
        done
        rm -f "$tmp_dir/$efile"
    done
    rm -rf "$tmp_dir" # Ensure temp dir is removed

    # --- CONSOLIDATION / RENAMING ---
    echo "[Member $mem] Finalizing..."
    # Extract unique variable/station pairs for this member
    targets=$(ls PROC_"${mem}"_*_st*.ts 2>/dev/null | sed "s/.*_\(.*\)_st\([0-9]*\)\.ts/\1_st\2/" | sort -u)

    for tg in $targets; do
        output="TOTAL_${mem}_${tg}.ts"
        if [ "$MODE" == "AN" ]; then
            # Merge multiple time fragments (an001, an002...)
            ls PROC_"${mem}"_an*_"${tg}".ts 2>/dev/null | sort -V | xargs cat > "$output"
        else
            # Single file, just rename
            mv PROC_"${mem}"_"${INPUT_NAME}_${mem}.ext"_"${tg}".ts "$output"
        fi
        # Cleanup any remaining PROC files for this target
        rm -f PROC_"${mem}"_*_"${tg}".ts
    done
}

# --- PARALLEL EXECUTION ---
count=0
for mem in $members; do
    process_member "$mem" &
    ((count++))
    if [ "$count" -ge "$MAX_JOBS" ]; then
        wait
        count=0
    fi
done
wait

# --- STATISTICS PHASE ---
echo "--- Calculating Global Ensemble Statistics ---"
all_targets=$(ls TOTAL_*_st*.ts 2>/dev/null | sed 's/TOTAL_[^_]*_\(.*\)\.ts/\1/' | sort -u)

for target in $all_targets; do
    member_files=$(ls TOTAL_*_"${target}".ts 2>/dev/null | sort -V)
    n_files=$(echo "$member_files" | wc -w)

    if [ "$n_files" -gt 0 ]; then
        output_stats="STATS_ENSEMBLE_${target}.ts"
        echo "Generating $output_stats ($n_files members)"
        paste $member_files | awk -v n="$n_files" '
        {
            sum=0; sumsq=0; max=-1e30; min=1e30;
            timestamp=$1;
            for (i=0; i<n; i++) {
                val=$( (i*2)+2 );
                sum+=val; sumsq+=val*val;
                if(val>max) max=val; if(val<min) min=val;
            }
            mean=sum/n;
            var=(sumsq/n)-(mean*mean);
            std=(var>0)?sqrt(var):0;
            printf "%s\t%.6f\t%.6f\t%.6f\t%.6f\n", timestamp, mean, std, max, min;
        }' > "$output_stats"
    fi
done

mkdir -p ext_results
mv TOTAL_*.ts STATS_ENSEMBLE_*.ts ext_results

echo "Workflow completed successfully."
