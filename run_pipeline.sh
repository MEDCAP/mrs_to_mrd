#!/usr/bin/env bash
#
# Convert, shift-correct, reconstruct and plot every experiment subfolder of one directory:
#
#   MRStomrd2.py -t <name>.tar -o <name>_raw.mrd2   one tar per subfolder, one stream out
#   mrd2shift.py --method contiguous             -> <name>_contiguous.mrd2
#   mrd2recon.py on the raw and the shifted file -> <name>_raw_recon.mrd2,
#                                                   <name>_contiguous_recon.mrd2
#   mrdplot.py -s <target>/png on both recons    -> the figures, and a .mat, per recon
#
#   ./run_pipeline.sh <target_directory>
#
# Conversion runs on a tar rather than on the folder because -t is the only mode that takes -o:
# `organize_members` groups one archive into one stream, so the output path is named here rather
# than by the converter, and every later stage is handed an exact file. Everything lands beside
# the subfolders in the target directory, and the tar is removed once it has been converted.
#
# A subfolder that fails any stage is reported and abandoned there, and the sweep carries on
# with the next one - hence no `set -e`, which would take the whole run down with it.

set -uo pipefail

# 1. Validate argument count
if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <target_directory>" >&2
  exit 1
fi

TARGET_DIR="$1"

# 2. Validate input is a directory
if [ ! -d "$TARGET_DIR" ]; then
  echo "Error: Directory '$TARGET_DIR' does not exist." >&2
  exit 1
fi

# Resolve to an absolute path before leaving for the script directory
TARGET_DIR="$(cd "$TARGET_DIR" && pwd)"
OUTPUT_DIR="${TARGET_DIR}/png"
mkdir -p "$OUTPUT_DIR"

# Ensure running in script directory so Python scripts are resolved
cd "$(cd "$(dirname "$0")" >/dev/null 2>&1 && pwd)"

export MPLBACKEND=Agg

# Helper function to convert MSYS/Git Bash paths for native Windows Python
winpath() {
  if command -v cygpath >/dev/null 2>&1; then
    cygpath -w "$1"
  else
    echo "$1"
  fi
}

# Reconstruct one stream and plot what came out. Takes the .mrd2 to reconstruct; the recon file
# is named after it, so the raw and the shifted recon sit side by side under their own names.
recon_and_plot() {
  local stream_path="$1"
  local recon_path="${stream_path%.mrd2}_recon.mrd2"

  echo "=== $(basename "$stream_path"): recon ==="
  if ! python mrd2recon.py -i "$(winpath "$stream_path")" -o "$(winpath "$recon_path")" \
       -bic_tm 0.0 -urea 2.3 -pyr_s 9.7 -ala_tm 15.2 -hyd_tm 18.1 -lac_m 21.8; then
    echo "Error: $(basename "$stream_path") failed to reconstruct." >&2
    rm -f "$recon_path"
    return 1
  fi

  echo "=== $(basename "$recon_path"): plotting ==="
  if ! python mrdplot.py -i "$(winpath "$recon_path")" -s "$(winpath "$OUTPUT_DIR")"; then
    echo "Error: $(basename "$recon_path") failed to plot." >&2
    return 1
  fi
  return 0
}

# 3. Iterate through immediate subdirectories
shopt -s nullglob
found_subdirs=0
failed=()

for dir in "$TARGET_DIR"/*/; do
  [ -d "$dir" ] || continue

  # Skip OUTPUT_DIR if located inside TARGET_DIR
  if [ "$(cd "$dir" && pwd)" = "$OUTPUT_DIR" ]; then
    continue
  fi

  found_subdirs=$((found_subdirs + 1))

  dir="${dir%/}"
  base_name="$(basename "$dir")"
  archive_path="${TARGET_DIR}/${base_name}.tar"
  raw_path="${TARGET_DIR}/${base_name}_raw.mrd2"
  shift_path="${TARGET_DIR}/${base_name}_contiguous.mrd2"

  # one root directory in the archive, which is what names the measurement inside the stream
  echo "=== ${base_name}: archiving ==="
  if ! tar cf "$archive_path" -C "$TARGET_DIR" "$base_name"; then
    echo "Error: ${base_name} could not be archived; moving on." >&2
    rm -f "$archive_path"
    failed+=("$base_name")
    continue
  fi

  echo "=== ${base_name}: converting ==="
  if ! python MRStomrd2.py -t "$(winpath "$archive_path")" -o "$(winpath "$raw_path")"; then
    echo "Error: ${base_name} failed to convert; moving on." >&2
    rm -f "$archive_path" "$raw_path"
    failed+=("$base_name")
    continue
  fi
  rm -f "$archive_path"

  echo "=== ${base_name}: shift correction ==="
  if ! python mrd2shift.py -i "$(winpath "$raw_path")" -o "$(winpath "$shift_path")" \
       --method contiguous; then
    echo "Error: ${base_name} failed the shift correction; moving on." >&2
    rm -f "$shift_path"
    failed+=("$base_name")
    continue
  fi

  # the uncorrected stream and the shifted one beside it, so the two can be compared
  subfolder_ok=1
  for stream_path in "$raw_path" "$shift_path"; do
    if ! recon_and_plot "$stream_path"; then
      subfolder_ok=0
      break
    fi
  done
  if [ "$subfolder_ok" -eq 0 ]; then
    echo "Error: ${base_name} did not finish; moving on." >&2
    failed+=("$base_name")
  fi
done

if [ "$found_subdirs" -eq 0 ]; then
  echo "No subdirectories found in '$TARGET_DIR'."
  exit 1
fi
if [ "${#failed[@]}" -ne 0 ]; then
  echo "Done. ${#failed[@]} of ${found_subdirs} subdirectories failed: ${failed[*]}."
  echo "Figures are in ${OUTPUT_DIR}."
  exit 1
fi
echo "Done. ${found_subdirs} subdirectories processed; figures are in ${OUTPUT_DIR}."
