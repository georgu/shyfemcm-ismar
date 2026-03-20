#!/usr/bin/env bash
#------------------------------------------------------------------------------
# SHY CONCATENATOR (Serial Version)
# Concatenates all an?????_enNNNNNb.[type].shy files for each member.
# Supports both 'hydro.shy' and 'ts.shy' separately.
# Result: enNNNNN.hydro.shy and enNNNNN.ts.shy
#------------------------------------------------------------------------------

# --- CONFIGURATION ---
SCRIPTPATH="$(cd "$(dirname "$0")" && pwd)"
FEMDIR="${SCRIPTPATH}/../../.."
SHYELAB="${FEMDIR}/bin/shyelab"

# Check dependencies
[ -x "$SHYELAB" ] || { echo "FATAL: shyelab not found at $SHYELAB"; exit 1; }

# Types of SHY files to look for
shy_types="hydro.shy ts.shy"

# --- MAIN LOGIC ---

# 1. Identify unique ensemble members (extracting NNNNN from any matching file)
shopt -s nullglob
candidates=( an*_en*.*.shy )
shopt -u nullglob

if [ ${#candidates[@]} -eq 0 ]; then
  echo "No .shy files found matching the pattern an*_en*.*.shy"
  exit 1
fi

# Build unique, sorted member list
members=$(ls an*_en*.*.shy 2>/dev/null | sed -n 's/.*_en\([0-9]*\)b\..*/\1/p' | sort -u)

# 2. Sequential Loop (Serial)
for nen in $members; do
  echo "=== Processing Member: $nen ==="

  for stype in $shy_types; do
    shopt -s nullglob
    # Find files like an00001_en00001b.hydro.shy OR an00001_en00001b.ts.shy
    files=( an*_en"${nen}"b."${stype}" )
    shopt -u nullglob

    if [ ${#files[@]} -gt 0 ]; then
      echo "Concatenating ${#files[@]} files of type [${stype}] for member $nen..."

      # shyelab concatenation (catmode +1)
      # Output to temporary out.shy
      "$SHYELAB" -out "out.shy" -catmode +1 "${files[@]}" > /dev/null 2>&1

      if [ -f "out.shy" ]; then
        # Final name: en00001.hydro.shy or en00001.ts.shy
        mv -f "out.shy" "en${nen}.${stype}"
        echo "Done: en${nen}.${stype} created."
      else
        echo "ERROR: shyelab failed to create out.shy for member $nen type $stype"
      fi
    else
      # If the type (e.g. ts.shy) is missing for this member, just skip it
      continue
    fi
  done
  echo "------------------------------------------------"
done

echo "All SHY merges (hydro and ts) completed."
