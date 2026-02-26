#!/usr/bin/env bash
#------------------------------------------------------------------------------
# Run many SHYFEM simulations using different .str files with a common basename.
# Uses GNU Parallel.
#
# Usage: shyfem_ens.sh [n_threads] [basename]
#   n_threads : integer (0 -> use as many jobs as .str files)
#   basename  : common prefix of .str files to run (matches "[basename]*.str")
#------------------------------------------------------------------------------

Usage() {
  echo
  echo "Run many SHYFEM simulations using different str files with a common"
  echo "basename. This program uses GNU Parallel."
  echo
  echo "Usage: $(basename "$0") [n_threads] [basename]"
  echo "  n_threads : 0 uses the number of .str files (cap to available jobs)"
  echo "  basename  : base name of the .str files (pattern: \"[basename]*.str\")"
  exit 0
}

die() { echo "FATAL: $*" >&2; exit 1; }

#------------------------------------------------------------------------------
# Locate script, SCRIPTPATH, FEMDIR (robust even without realpath)
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
FEMDIR="${SCRIPTPATH}/.."   # fem directory (as in the original)

#------------------------------------------------------------------------------
# Args parsing and validation
#------------------------------------------------------------------------------
[ $# -eq 2 ] || Usage
nth="$1"
fbasename="$2"

# integer check for nth
[[ "$nth" =~ ^-?[0-9]+$ ]] || die "n_threads must be an integer, got '$nth'."

#------------------------------------------------------------------------------
# Dependency checks
#------------------------------------------------------------------------------
command -v parallel >/dev/null 2>&1 || die "'parallel' is not installed."
[ -x "${FEMDIR}/fem3d/shyfem" ] || die "Executable not found: ${FEMDIR}/fem3d/shyfem"

#------------------------------------------------------------------------------
# Collect .str files matching the basename (safe glob, no 'ls' parsing)
#------------------------------------------------------------------------------
shopt -s nullglob
strfiles=( "${fbasename}"*.str )
shopt -u nullglob
[ "${#strfiles[@]}" -gt 0 ] || die "No .str files found matching '${fbasename}*.str'."

# basic existence/size check for each .str
for strf in "${strfiles[@]}"; do
  [ -s "$strf" ] || die "File '$strf' does not exist or has zero size."
done

#------------------------------------------------------------------------------
# Decide parallelism: if nth<=0, use number of files; cap to number of files
#------------------------------------------------------------------------------
if [ "$nth" -le 0 ] || [ "$nth" -gt "${#strfiles[@]}" ]; then
  nth="${#strfiles[@]}"
fi

#------------------------------------------------------------------------------
# Worker: run a single simulation
#------------------------------------------------------------------------------
Make_sim() {
  local strfile="$1" fem3d_dir="$2"
  local basen
  basen="$(basename "$strfile" .str)"
  "${fem3d_dir}/shyfem" "$strfile" > "${basen}.log"
}
export -f Make_sim
export FEMDIR

#------------------------------------------------------------------------------
# Run all jobs in parallel
#------------------------------------------------------------------------------
echo "Running ${#strfiles[@]} simulations with parallelism P=${nth} ..."
# shellcheck disable=SC2086
parallel --no-notice -P "$nth" Make_sim ::: "${strfiles[@]}" ::: "${FEMDIR}/fem3d"

echo "All simulations completed."
exit 0
