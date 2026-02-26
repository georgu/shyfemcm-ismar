#!/usr/bin/env bash
#------------------------------------------------------------------------------
# Copyright (C) 2017
# Marco Bajo, CNR-ISMAR Venice. All rights reserved.
#
# Split SHYFEM .ext files into per-variable fragments and rename them as .ts.
# Usage: split_exts.sh [basename]
# Processes all files matching '[basename]*.ext'
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
FEMDIR="${SCRIPTPATH}/.."   # fem directory (as in the original)

Usage() {
  echo "Usage: $(basename "$0") [basename]"
  echo
  echo "Splits all files matching '[basename]*.ext' with shyelab -split and"
  echo "renames parts into '<base>_<var.dim.level>.ts'."
  exit 0
}

die() { echo "FATAL: $*" >&2; exit 1; }

#------------------------------------------------------------------------------
# Args
#------------------------------------------------------------------------------
[ $# -eq 1 ] || Usage
fbase="$1"

#------------------------------------------------------------------------------
# Dependencies
#------------------------------------------------------------------------------
[ -x "${FEMDIR}/fembin/shyelab" ] || die "shyelab not found at ${FEMDIR}/fembin/shyelab."

#------------------------------------------------------------------------------
# Variables list to look for after split (kept from your original script)
#------------------------------------------------------------------------------
vars='zeta.2d velx.2d vely.2d velx.3d vely.3d speed.2d speed.3d dir.2d dir.3d all.2d temp.2d temp.3d salt.2d salt.3d'

#------------------------------------------------------------------------------
# Expand target files safely (no 'ls' parsing)
#------------------------------------------------------------------------------
shopt -s nullglob
ext_files=( "${fbase}"*.ext )
shopt -u nullglob

[ "${#ext_files[@]}" -gt 0 ] || die "No '.ext' files found matching prefix '${fbase}'."

#------------------------------------------------------------------------------
# Process each .ext file
#------------------------------------------------------------------------------
for efile in "${ext_files[@]}"; do
  echo "File: $efile"
  basefile="$(basename "$efile" .ext)"

  # Split into var.* chunks in the current directory
  # Log to base-specific log to avoid clobbering
  "${FEMDIR}/fembin/shyelab" -split "$efile" > "${basefile}.ext.log" 2>&1

  # Remove previous outputs for this base (only those produced by this step)
  rm -f "${basefile}_"*.ts || true

  # For each variable group, rename all generated levels var.* -> <base>_var.*.ts
  for vv in $vars; do
    # If there is at least one level for this variable, process them
    shopt -s nullglob
    vfiles=( ${vv}.* )
    shopt -u nullglob

    [ "${#vfiles[@]}" -gt 0 ] || continue

    for fl in "${vfiles[@]}"; do
      mv -f "$fl" "${basefile}_${fl}.ts"
    done
  done
done

echo "Done."
