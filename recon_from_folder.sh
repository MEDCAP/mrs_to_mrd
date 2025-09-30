#!/bin/bash

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
echo "run reconstruction code"
winpty python mrd2recon.py -f "$folder_path" -p "-bic 0.0 -urea 2.3 -pyr 9.7 -ala 15.2 -poop 15.9 -hyd 18.1 -lac 21.8"
# winpty python mrdplot.py -i "$folder_path\recon.mrd2"