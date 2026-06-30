#!/bin/bash

# Exit on error and trap errors to keep window open
set -e
trap 'echo ""; echo "❌ Script failed! Press any key to exit..."; read -n 1' ERR

# ---------------------------------------------------------------------------
# 1. Determine the data folder (cross-platform)
#    Priority: CLI argument > Windows PowerShell folder dialog
# ---------------------------------------------------------------------------
folder_path=""
if [ -n "$1" ]; then
    folder_path="$1"
elif command -v powershell.exe &> /dev/null; then
    # A single-line PowerShell command to open the dialog.
    # Note the backslashes to escape the dollar signs for bash.
    win_path=$(powershell.exe -Command "
    Add-Type -AssemblyName System.Windows.Forms;
    \$folderBrowser = New-Object System.Windows.Forms.FolderBrowserDialog;
    \$folderBrowser.Description = 'Select a folder';
    if (\$folderBrowser.ShowDialog() -eq 'OK') {
        \$folderBrowser.SelectedPath
    }")
    if [ -n "$win_path" ]; then
        # If using WSL, convert the path from C:\... to /mnt/c/...
        if command -v wslpath &> /dev/null; then
            folder_path=$(wslpath "$win_path")
        else
            # For Git Bash, just use the Windows path.
            folder_path="$win_path"
        fi
    fi
fi

if [ -z "$folder_path" ]; then
    echo "❌ No folder was selected. Pass a folder path as the first argument."
    echo "Press any key to exit..."
    read -n 1
    exit 1
fi
echo "You selected: $folder_path"

# ---------------------------------------------------------------------------
# 2. Cross-platform python runner (winpty is only present on Windows/Git Bash)
# ---------------------------------------------------------------------------
if command -v winpty &> /dev/null; then
    PY="winpty python"
else
    PY="python"
fi

# ---------------------------------------------------------------------------
# 3. Output directories. On macOS/Linux point the recon save dirs at the local
#    data tree; on Windows leave them unset so mrd2recon.py uses its defaults
#    (C:/Users/MRS/Desktop/shurik/...).
# ---------------------------------------------------------------------------
case "$(uname -s)" in
    Darwin|Linux)
        export MRS_PROCESSED_DIR="$HOME/dev/data/epsi_kidney_data/processed_npyfiles"
        export MRS_LORNFIT_DIR="$HOME/dev/data/epsi_kidney_data/lorn_fit"
        ;;
    *)
        : # Windows: rely on mrd2recon.py defaults
        ;;
esac
mat_dir="${MRS_PROCESSED_DIR:-C:/Users/MRS/Desktop/shurik/processed_npyfiles}"
mkdir -p "$mat_dir" 2>/dev/null || true
if [ -n "$MRS_LORNFIT_DIR" ]; then
    mkdir -p "$MRS_LORNFIT_DIR" 2>/dev/null || true
fi

# ---------------------------------------------------------------------------
# 4. Convert MRS to MRD format (once). Already-converted folders are skipped
#    by MRStomrd2.py (they already contain raw.mrd2).
# ---------------------------------------------------------------------------
echo "Converting MRS to MRD format..."
if ! $PY MRStomrd2.py -f "$folder_path" -u 3; then
    echo "❌ MRS to MRD conversion failed!"
    echo "Press any key to exit..."
    read -n 1
    exit 1
fi

# ---------------------------------------------------------------------------
# 5. Read the list of newly-converted experiments from the manifest that
#    MRStomrd2.py writes (one raw.mrd2 path per line). Only these folders are
#    reconstructed; already-converted data is skipped.
#    (while-read keeps this compatible with bash 3.2, the macOS default.)
# ---------------------------------------------------------------------------
manifest_file="$folder_path/newly_converted_mrd_files.txt"
if [ ! -f "$manifest_file" ]; then
    echo "❌ No convert manifest found at $manifest_file"
    echo "Press any key to exit..."
    read -n 1
    exit 1
fi

raw_files=()
while IFS= read -r raw; do
    # skip blank lines
    [ -n "${raw//[[:space:]]/}" ] && raw_files+=("$raw")
done < "$manifest_file"

if [ ${#raw_files[@]} -eq 0 ]; then
    echo "No newly converted experiments listed in $manifest_file"
    echo "Press any key to exit..."
    read -n 1
    exit 1
fi

echo "Found ${#raw_files[@]} newly converted experiment(s) to reconstruct."

# ---------------------------------------------------------------------------
# 6. Reconstruct each experiment twice:
#      RUN1: with poop  (7 peaks) -> <exp>_metaboites_withpoop.mat
#      RUN2: without poop (6 peaks) -> <exp>_metabolites.mat
#    Delete *_recon.mrd2 before each pass so the second run is not skipped
#    by mrd2recon.py's "recon already exists" check.
# ---------------------------------------------------------------------------
RUN1=(-bic_tm 0.0 -urea 2.3 -pyr_s 9.7 -ala_tm 15.2 -poop_tm 15.9 -hyd_tm 18.1 -lac_m 21.8)
RUN2=(-bic_tm 0.0 -urea 2.3 -pyr_s 9.7 -ala_tm 15.2 -hyd_tm 18.1 -lac_m 21.8)

for raw in "${raw_files[@]}"; do
    exp_dir=$(dirname "$raw")
    exp=$(basename "$exp_dir")

    echo "=== $exp : run 1 (with poop) ==="
    rm -f "$exp_dir"/*_recon.mrd2
    $PY mrd2recon.py -f "$raw" "${RUN1[@]}" || echo "⚠️  Reconstruction (run 1) failed for $exp"

    echo "=== $exp : run 2 (without poop) ==="
    rm -f "$exp_dir"/*_recon.mrd2
    $PY mrd2recon.py -f "$raw" "${RUN2[@]}" || echo "⚠️  Reconstruction (run 2) failed for $exp"
done

# ---------------------------------------------------------------------------
# 7. Confirm both mat files exist for each experiment and write a report into
#    the mat output folder.
# ---------------------------------------------------------------------------
report="$mat_dir/matfile_check.txt"
: > "$report"
all_ok=1

for raw in "${raw_files[@]}"; do
    exp=$(basename "$(dirname "$raw")")
    withpoop="$mat_dir/${exp}_metaboites_withpoop.mat"
    plain="$mat_dir/${exp}_metabolites.mat"
    if [ -f "$withpoop" ] && [ -f "$plain" ]; then
        echo "$exp: OK (both mat files present)" >> "$report"
    else
        all_ok=0
        missing=""
        [ -f "$withpoop" ] || missing="$missing _metaboites_withpoop.mat"
        [ -f "$plain" ]    || missing="$missing _metabolites.mat"
        echo "$exp: MISSING -$missing" >> "$report"
    fi
done

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Mat file check report written to: $report"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
cat "$report"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

if [ "$all_ok" -eq 1 ]; then
    echo "✅ All experiments have both mat files."
else
    echo "⚠️  Some mat files are missing (see report above)."
fi

# ---------------------------------------------------------------------------
# 8. Clean up the convert manifest text file generated during this run
#    (the matfile_check.txt report in the mat folder is kept as the deliverable).
# ---------------------------------------------------------------------------
rm -f "$manifest_file"
echo "Cleaned up manifest: $manifest_file"
