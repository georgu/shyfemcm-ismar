#!/usr/bin/env bash
#------------------------------------------------------------------------------
# Copyright (C) 2017
# Marco Bajo, CNR-ISMAR Venice. All rights reserved.
#
# Merge of ensemble outputs over an analysis period.
# - file_type = 'ext': split by variables and merge to timeseries .ts per var/level
# - file_type = 'shy': concatenate .shy hydrodynamic files into a single shy
#
# Usage: merge_ens.sh [file_type] [output]
#   file_type : ext | shy
#   output    : small (z) | medium (z,t,s) | full (all)   # only for 'ext'
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Resolve script path (fallback if realpath is missing)
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
FEMDIR="${SCRIPTPATH}/.."   # fem directory
SIMDIR="$(pwd)"             # current dir

Usage() {
  cat <<EOF

Usage: $(basename "$0") [file_type] [output]

  file_type : ext | shy
  output    : small | medium | full      (ONLY for 'ext')

Notes:
  - For 'ext': merges split variables into timeseries .ts files:
       small  -> zeta.2d
       medium -> zeta.2d temp.{2d,3d} salt.{2d,3d}
       full   -> a curated set of 2d/3d fields (see code)
  - For 'shy': concatenates all an?????_enNNNNNb.hydro.shy for each member.

EOF
  exit 0
}

die() { echo "FATAL: $*" >&2; exit 1; }
require_cmd() { command -v "$1" >/dev/null 2>&1 || die "Required command '$1' not found."; }
check_file() { [ -s "$1" ] || die "File '$1' does not exist or has zero size."; }

#------------------------------------------------------------------------------
# Dependency checks
#------------------------------------------------------------------------------
require_cmd awk
require_cmd sed
require_cmd cut
[ -x "${FEMDIR}/fembin/shyelab" ] || die "shyelab not found at ${FEMDIR}/fembin/shyelab."

#------------------------------------------------------------------------------
# Parse args
#------------------------------------------------------------------------------
if [ $# -ge 1 ]; then
  file_type="$1"
else
  Usage
fi

outt="${2:-}"   # only used for ext

case "$file_type" in
  ext)
    case "$outt" in
      small|medium|full) : ;;
      *) Usage ;;
    esac
    ;;
  shy)
    : # 'outt' ignored
    ;;
  *)
    Usage
    ;;
esac

#------------------------------------------------------------------------------
# Merge EXT: produce per-variable, per-level .ts files for each ensemble member
#------------------------------------------------------------------------------
Merge_timeseries() {
  local nen="$1" ftype="$2" outt="$3"

  # Collect all ext files matching this member (e.g., an*_en00012b.ext)
  shopt -s nullglob
  local files=( an*_en"${nen}"b."${ftype}" )
  shopt -u nullglob
  [ "${#files[@]}" -gt 0 ] || die "No '${ftype}' files found for member ${nen}."

  # Variable lists (kept compatible with your original)
  local allvars="all.2d dir.2d salt.2d speed.2d temp.2d velx.2d vely.2d zeta.2d dir.3d salt.3d speed.3d temp.3d velx.3d vely.3d"
  local vars
  case "$outt" in
    small)  vars="zeta.2d" ;;
    medium) vars="zeta.2d temp.2d salt.2d temp.3d salt.3d" ;;
    full)   vars="$allvars" ;;
    *) die "Invalid output selector '$outt'." ;;
  esac

  # Clean previous partials for this member
  rm -f *_st*_en"${nen}".ts *.2d.* *.3d.* || true

  # Process each ext file: split and append selected variables
  for fil in "${files[@]}"; do
    check_file "$fil"
    echo "Processing file: $fil"
    # Split into <var>.<lev>.ext / shyelab writes var.* in cwd
    "${FEMDIR}/fembin/shyelab" -split "$fil" &>/dev/null

    # For each selected variable, scan all its levels and append (drop header/footer)
    for vv in $vars; do
      shopt -s nullglob
      # files like temp.3d.0001 (var.dim.level)
      for flev in ${vv}.*; do
        # Extract level index: third field after dots
        local idst
        idst="$(echo "$flev" | cut -d '.' -f 3)"
        # Append content without first and last line
        # (equiv. to your intent, but pipeline fixed and robust)
        if [ -s "$flev" ]; then
          # ensure output file name (consistent with your naming)
          local out_ts="${vv}_st${idst}_en${nen}.ts"
          # drop first header line and last footer line if present
          # guard for 2-line files: will yield empty (which is fine)
          tail -n +2 "$flev" | head -n -1 >> "$out_ts" || true
        fi
        rm -f "$flev"
      done
      shopt -u nullglob
    done
  done

  # final cleanup of any residual split files
  rm -f *.2d.* *.3d.* || true
}

#------------------------------------------------------------------------------
# Merge SHY: concatenate hydrodynamic shy files for a member
#------------------------------------------------------------------------------
Merge_shy() {
  local nen="$1"

  shopt -s nullglob
  local files=( an*_en"${nen}"b.hydro.shy )
  shopt -u nullglob
  [ "${#files[@]}" -gt 0 ] || die "No '.hydro.shy' files found for member ${nen}."

  # shyelab concatenation (catmode +1)
  "${FEMDIR}/fembin/shyelab" -out -catmode +1 "${files[@]}"
  mv -f out.shy "en${nen}.hydro.shy"
}

#------------------------------------------------------------------------------
# Enumerate ensemble members by scanning any step with a stable pattern
# Your original used: `for efile in $(ls an00002_en*.${file_type})` and substring
# Now: robust glob + parsing of the member id between 'en' and 'b'
#------------------------------------------------------------------------------
shopt -s nullglob
candidates=( an*_en*."${file_type}" )
shopt -u nullglob
[ "${#candidates[@]}" -gt 0 ] || die "No input files found for type '${file_type}'."

# Build a unique, sorted list of member IDs (zero-padded 5 digits)
declare -A seen
members=()
for efile in "${candidates[@]}"; do
  # match pattern an?????_enNNNNNb.<ext or shy...>
  bn="$(basename "$efile")"
  # Extract part after 'en' and before 'b.'
  # Example: an00002_en00012b.ext -> id=00012
  id="$(echo "$bn" | sed -n 's/.*_en\([0-9][0-9][0-9][0-9][0-9]\)b\..*/\1/p')"
  [ -n "${id:-}" ] || die "Cannot parse member id from '$bn'."
  if [ -z "${seen[$id]:-}" ]; then
    seen[$id]=1
    members+=( "$id" )
  fi
done
IFS=$'\n' members=( $(printf "%s\n" "${members[@]}" | sort) ); unset IFS

#------------------------------------------------------------------------------
# Loop over members and merge according to file_type
#------------------------------------------------------------------------------
for nen in "${members[@]}"; do
  echo "=== Merging for ensemble member ${nen} (${file_type}) ==="
  case "$file_type" in
    ext) Merge_timeseries "$nen" "ext" "$outt" ;;
    shy) Merge_shy        "$nen" ;;
  esac
done

echo "All merges completed."
exit 0
