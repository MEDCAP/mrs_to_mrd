#!/bin/bash

# 1. Select folder to feed to recon with slicing arguments or select
#    If parent folder is passed, 
# 2. Ideally use the 



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

    echo "✅ You selected: $folder_path"
    echo "--- Contents of the folder ---"
    ls -la "$folder_path"
else
    echo "❌ No folder was selected."
fi

# call python scripts
winpty python MRStomrd2.py -f "$folder_path" -u 3 
echo "Conversion to MRD2 completed."
winpty python mrd2recon.py -f "$folder_path" -bic_tm 0.0 -urea 2.3 -pyr_s 9.7 -ala_tm 15.2 -poop_tm 15.9 -hyd_tm 18.1 -lac_m 21.8

echo "🔍 Searching for recon.mrd2 files..."

# Find all files that end with recon.mrd2 files under the folder path
mapfile -t mrd2_files < <(find "$folder_path" -name "*recon.mrd2" -type f 2>/dev/null)

# Check if any files were found
if [ ${#mrd2_files[@]} -eq 0 ]; then
    echo "❌ No recon.mrd2 files found in $folder_path"
    exit 1
fi

echo "📁 Found ${#mrd2_files[@]} recon.mrd2 file(s):"
echo ""

# Display numbered list of files
for i in "${!mrd2_files[@]}"; do
    echo "$((i+1)). ${mrd2_files[i]}"
done

echo ""
echo "Please select a file by entering its number:"
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