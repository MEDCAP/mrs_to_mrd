#!/bin/bash

# Exit on error and trap errors to keep window open
set -e
trap 'echo ""; echo "❌ Script failed! Press any key to exit..."; read -n 1' ERR

# A single-line PowerShell command to open the dialog.
# Note the backslashes to escape the dollar signs for bash.
win_path=$(powershell.exe -Command "
Add-Type -AssemblyName System.Windows.Forms;
\$folderBrowser = New-Object System.Windows.Forms.FolderBrowserDialog;
\$folderBrowser.Description = 'Select a folder';
if (\$folderBrowser.ShowDialog() -eq 'OK') {
    \$folderBrowser.SelectedPath
}")

# Check if a path was returned.
if [ -n "$win_path" ]; then
    # If using WSL, convert the path from C:\... to /mnt/c/...  
    if command -v wslpath &> /dev/null; then
        folder_path=$(wslpath "$win_path")
    else
        # For Git Bash, you might need a different conversion or just use the Windows path.
        folder_path="$win_path"
    fi
    echo "You selected: $folder_path"

else
    echo "❌ No folder was selected."
    echo "Press any key to exit..."
    read -n 1
    exit 1
fi

# call python scripts
echo "Converting MRS to MRD format..."
if ! winpty python MRStomrd2.py -f "$folder_path" -u 1; then
    echo "❌ MRS to MRD conversion failed!"
    echo "Press any key to exit..."
    read -n 1
    exit 1
fi

echo "Running reconstruction code..."
if ! winpty python mrd2recon.py -f "$folder_path" -urea 0.0 -KIC_s 8.6 -leu_tm 13.0 -hyd_tm 18.1 -?_tm 21.8 -w 1.; then
    echo "❌ Reconstruction failed!"
    echo "Press any key to exit..."
    read -n 1
    exit 1
fi

echo "🔍 Searching for recon.mrd2 files..."

# Find all recon.mrd2 files under the folder path
mapfile -t mrd2_files < <(find "$folder_path" -name "*recon.mrd2" -type f 2>/dev/null)

# Check if any files were found
if [ ${#mrd2_files[@]} -eq 0 ]; then
    echo "No recon.mrd2 files found in $folder_path"
    echo "Press any key to exit..."
    read -n 1
    exit 1
fi

echo "Found ${#mrd2_files[@]} recon.mrd2 file(s):"
echo ""

# Add Ctrl+C handler for graceful exit
trap 'echo ""; echo "Exiting..."; exit 0' INT

# Loop until user presses Ctrl+C
while true; do
    echo ""
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "Available recon.mrd2 files:"
    echo ""
    
    # Display numbered list of files
    for i in "${!mrd2_files[@]}"; do
        echo "  $((i+1)). ${mrd2_files[i]}"
    done
    
    echo ""
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "Press Ctrl+C to exit"
    read -p "Select file number (1-${#mrd2_files[@]}): " selection
    
    # Validate input
    if ! [[ "$selection" =~ ^[0-9]+$ ]] || [ "$selection" -lt 1 ] || [ "$selection" -gt ${#mrd2_files[@]} ]; then
        echo "❌ Invalid selection. Please enter a number between 1 and ${#mrd2_files[@]}."
        continue
    fi
    
    # Get selected file path
    selected_file="${mrd2_files[$((selection-1))]}"
    echo ""
    echo "✅ Selected: $selected_file"
    
    # Plot the selected file
    if ! winpty python mrdplot.py -i "$selected_file"; then
        echo "❌ Plotting failed!"
        echo "Press any key to continue..."
        read -n 1
        continue
    fi
    
    echo ""
    echo "✅ Plotting complete!"
done