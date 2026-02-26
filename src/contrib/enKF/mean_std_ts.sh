#!/usr/bin/env bash
#------------------------------------------------------------------------------
# Copyright (C) 2021
# Marco Bajo, CNR-ISMAR Venice. All rights reserved.
#
# Compute mean and std of an ensemble of time series (columns: time value).
# All files must have the same number of records; timestamps are assumed aligned.
# Output: <basename>mean_std.dat  (time, mean, std)
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Resolve script path (with fallback if realpath is missing)
#------------------------------------------------------------------------------
if command -v realpath >/dev/null 2>&1; then
  SCRIPT="$(realpath "$0")"
else
  SCRIPT="$0"
  while [ -L "$SCRIPT" ]; do
    LINK="$(readlink "$SCRIPT")"
    case "$LINK" in
      /*) SCRIPT="$LINK" ;;
      *)  SCRIPT="$(cd "$(dirname "$SCRIPT")" && pwd)/$LINK" ;;
    esac
  done
  SCRIPT="$(cd "$(dirname "$SCRIPT")" && pwd)/$(basename "$SCRIPT")"
fi
SCRIPTPATH="$(cd "$(dirname "$SCRIPT")" && pwd)"
FEMDIR="${SCRIPTPATH}/.."     # fem directory (kept for compatibility)   # not used directly
SIMDIR="$(pwd)"               # current dir

Usage() {
  echo "Usage: ./mean_std_ts.sh [basename]"
  echo
  echo "Reads all files matching '[basename]*' and writes '[basename]mean_std.dat'."
  exit 0
}

#------------------------------------------------------------------------------
# Check prerequisites
#------------------------------------------------------------------------------
command -v awk >/dev/null 2>&1 || { echo "FATAL: 'awk' is required." >&2; exit 1; }

#------------------------------------------------------------------------------
# Args
#------------------------------------------------------------------------------
[ "$#" -eq 1 ] || Usage
bfile="$1"

# Expand the set of input files safely (no 'ls' parsing)
shopt -s nullglob
files=( "${bfile}"* )
shopt -u nullglob

# Validate file list
if [ "${#files[@]}" -eq 0 ]; then
  echo "FATAL: No files found matching prefix '${bfile}'." >&2
  exit 1
fi

# Check each file
k=0
for tfile in "${files[@]}"; do
  if [ ! -s "$tfile" ]; then
    echo "FATAL: file '$tfile' does not exist or has zero size." >&2
    exit 1
  fi
  echo "Ensemble member: $tfile"
  k=$((k+1))
done

echo
echo "Number of ensemble members: $k"
if [ "$k" -lt 3 ]; then
  echo "Stopping: need at least 3 members to compute a stable std." >&2
  exit 1
fi

# Ensure all files have the same number of rows (no check on time values)
# We take the first file as reference.
ref_rows="$(wc -l < "${files[0]}")"
for tfile in "${files[@]:1}"; do
  nrows="$(wc -l < "$tfile")"
  if [ "$nrows" -ne "$ref_rows" ]; then
    echo "FATAL: Inconsistent length: '${files[0]}' has ${ref_rows} rows, '$tfile' has ${nrows}." >&2
    exit 1
  fi
done

# Output file (same naming as the original script)
out="${bfile}mean_std.dat"

# Compute mean and std:
#  - T[j] = time of row j (from any file; lengths are equal)
#  - v1[j] = sum of values; v2[j] = sum of squares
#  - std = sqrt( max( v2/n - (v1/n)^2 , 0 ) + eps )
LC_ALL=C awk -v nr="$k" '
  {
    T[FNR] = $1
    v1[FNR] += $2
    v2[FNR] += ($2 * $2)
    if (FNR > jmax) jmax = FNR
  }
  END {
    eps = 1.0e-12
    for (j = 1; j <= jmax; j++) {
      m  = v1[j] / nr
      vv = v2[j] / nr - m * m
      if (vv < 0) vv = 0   # numerical guard
      s  = sqrt(vv + eps)
      printf("%s %.12g %.12g\n", T[j], m, s)
    }
  }
' "${files[@]}" > "$out"

echo "Wrote: $out"
