#!/bin/bash

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
else
    echo "❌ No folder was selected."
fi

# Find all files that end with recon.mrd2 files under the folder path
mapfile -t mrd2_files < <(find "$folder_path" -name "*recon.mrd2" -type f 2>/dev/null)

# Check if any files were found
if [ ${#mrd2_files[@]} -eq 0 ]; then
    echo "❌ No recon.mrd2 files found in $folder_path"
    exit 1
fi

echo "Found ${#mrd2_files[@]} recon.mrd2 file(s):"

while true; do
    # Display numbered list of files
    for i in "${!mrd2_files[@]}"; do
        echo "$((i+1)). ${mrd2_files[i]}"
    done

    echo ""
    echo "Please select a file by entering its number or Ctrl+C to exit:"
    read -p "Enter selection (1-${#mrd2_files[@]}): " selection

    # Validate input
    if ! [[ "$selection" =~ ^[0-9]+$ ]] || [ "$selection" -lt 1 ] || [ "$selection" -gt ${#mrd2_files[@]} ]; then
        echo "❌ Invalid selection. Please enter a number between 1 and ${#mrd2_files[@]}."
        exit 1
    fi

    # Get selected file path
    selected_file="${mrd2_files[$((selection-1))]}"
    echo ""
    echo "✅ Selected: $selected_file"

    # Use the selected file (example with your commented command)
    winpty python mrdplot.py -i "$selected_file"
done