#!/bin/bash

# Ensure both required arguments are provided
if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <pe_parameters_file> <base_template_file>"
  exit 1
fi

pe_file=$1
base_skel=$2

# Verify that the PE parameters file exists
if [[ ! -f "$pe_file" ]]; then
  echo "Error: PE parameters file '$pe_file' does not exist."
  exit 1
fi

# Verify that the base template file exists
if [[ ! -f "$base_skel" ]]; then
  echo "Error: Base template file '$base_skel' does not exist."
  exit 1
fi

echo "--- Starting flag restoration process ---"
echo "Using PE file: $pe_file"
echo "Using base template: $base_skel"

# 1. Read the PE parameters file line by line to extract the flags (e.g., 001_DCB)
while read -r flag_str dummy; do
  # Skip empty lines or lines starting with a comment character (#)
  [[ -z "$flag_str" || "$flag_str" == \#* ]] && continue

  # Find the exact line numbers where the flag appears in the base template file
  line_numbers=$(grep -n -F "$flag_str" "$base_skel" | cut -d':' -f1)

  # If the flag is not found in the base template, skip to the next one
  if [[ -z "$line_numbers" ]]; then
    echo "Warning: Flag '$flag_str' not found in $base_skel. Skipping."
    continue
  fi

  # Loop through each line number found (handles multiple occurrences of the same flag)
  for line in $line_numbers; do
    # Capture the exact content of the correct line from the base template
    correct_line=$(sed -n "${line}p" "$base_skel")

    echo "Found flag '$flag_str' at line $line. Restoring in member files..."

    # 2. Loop through all member files matching the pattern member_*.skel
    for skel_file in member_*.skel; do
      # Strictly skip the base template file and ensure the target file exists
      [[ "$skel_file" == "$base_skel" ]] && continue
      [[ ! -f "$skel_file" ]] && continue

      # Create a secure temporary file for safe editing
      tmp_file=$(mktemp)
      
      # Overwrite the specific line based on the line number index (NR)
      awk -v ln="$line" -v new_text="$correct_line" 'NR==ln {$0=new_text} {print}' "$skel_file" > "$tmp_file"
      
      # Replace the member file with the corrected version containing the text flags
      mv "$tmp_file" "$skel_file"
    done
  done

done < "$pe_file"

echo "--------------------------------------------------------"
echo "Process completed! All member_*.skel files have been updated with text flags for DA."

