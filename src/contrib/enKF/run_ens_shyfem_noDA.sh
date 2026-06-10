#!/usr/bin/env bash
#------------------------------------------------------------------------------
# Run SHYFEM simulations in parallel using GNU Parallel.
# Each simulation is based on a unique .str file.
#
# Usage: ./shyfem_ens.sh [n_parallel] [basename]
#   n_parallel : number of concurrent jobs (0 = all files in parallel)
#   basename   : prefix for .str files (e.g., 'run' matches run01.str, run02.str)
#------------------------------------------------------------------------------

# --- Functions ---
Usage() {
  cat <<EOF

Usage: $(basename "$0") [n_parallel] [basename]

  n_parallel : Number of concurrent simulations.
               Use 0 to run all matches at once (limited by CPU cores).
  basename   : Common prefix of the .str files (pattern: "[basename]*.str").

Example:
  $(basename "$0") 4 sensitivity_run
EOF
  exit 0
}

die() { echo "FATAL: $*" >&2; exit 1; }

# --- Path Resolution ---
if command -v realpath >/dev/null 2>&1; then
  SCRIPT="$(realpath "$0")"
else
  # Fallback for systems without realpath
  TARGET_FILE="$0"
  cd "$(dirname "$TARGET_FILE")" || exit 1
  TARGET_FILE="$(basename "$TARGET_FILE")"
  while [ -L "$TARGET_FILE" ]; do
    TARGET_FILE="$(readlink "$TARGET_FILE")"
    cd "$(dirname "$TARGET_FILE")" || exit 1
    TARGET_FILE="$(basename "$TARGET_FILE")"
  done
  SCRIPT="$(pwd -P)/$TARGET_FILE"
fi

SCRIPTPATH="$(cd "$(dirname "$SCRIPT")" && pwd)"
FEMDIR="${SCRIPTPATH}/../../.." # Adjust based on your directory structure
SHYFEM_BIN="${FEMDIR}/bin/shyfem"

# --- Argument Validation ---
[ $# -eq 2 ] || Usage

nth="$1"
fbasename="$2"

# Check if nth is a valid integer
[[ "$nth" =~ ^[0-9]+$ ]] || die "n_parallel must be a positive integer, got '$nth'."

# --- Dependency & File Checks ---
command -v parallel >/dev/null 2>&1 || die "GNU 'parallel' is required but not installed."
[ -x "$SHYFEM_BIN" ] || die "Executable not found or not executable: $SHYFEM_BIN"

# Collect .str files safely
shopt -s nullglob
strfiles=( "${fbasename}"*.str )
shopt -u nullglob

[ "${#strfiles[@]}" -gt 0 ] || die "No .str files found matching '${fbasename}*.str'."

# --- Decide Parallelism ---
# If nth is 0 or larger than file count, use file count
if [ "$nth" -le 0 ] || [ "$nth" -gt "${#strfiles[@]}" ]; then
  nth="${#strfiles[@]}"
fi

# --- Worker Function ---
# Exported so GNU Parallel can see it
run_shy_sim() {
  local strfile="$1"
  local bin_exe="$2"
  local log_name="${strfile%.str}.log"
  local err_name="${strfile%.str}.err"

  echo "Starting: $strfile"
  # Run shyfem and redirect output to log
  "$bin_exe" "$strfile" > "$log_name" 2>&1

  if [ $? -eq 0 ]; then
    echo "Completed: $strfile"
  else
    echo "FAILED: $strfile (Check $log_name)"
    [[ -e "fort.999" ]] && mv fort.999 $err_name
  fi
}
export -f run_shy_sim

# --- Execution ---
echo "----------------------------------------------------------------"
echo "Found ${#strfiles[@]} simulations matching prefix: $fbasename"
echo "Running with parallelism P = $nth"
echo "----------------------------------------------------------------"

# Run parallel
# --no-notice: suppresses the citation notice
# --jobs/-P: number of concurrent processes
parallel --no-notice --jobs "$nth" run_shy_sim {} "$SHYFEM_BIN" ::: "${strfiles[@]}"

echo "----------------------------------------------------------------"
echo "All parallel simulations completed."
