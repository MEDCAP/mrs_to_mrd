#!/bin/bash
# Read parent directory of the MRS data interactively
# Loop through all the files from 


read -p 

file_path = "$1"

# Function to get file creation date
get_creation_date() {
  local file_path="$1"

  if [ -f "$file_path" ]; then
    # Use stat to get the birth time (%w)
    # The output format can be customized
    creation_date=$(stat --format=%w "$file_path")

    if [ -n "$creation_date" ] && [ "$creation_date" != "-" ]; then
      echo "Creation date of '$file_path': $creation_date"
    else
      echo "Creation date not available for '$file_path' (filesystem may not support it)."
    fi
  else
    echo "Error: File '$file_path' not found."
  fi
}

# Example usage:
get_creation_date "$file_path"