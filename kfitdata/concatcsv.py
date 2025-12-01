import os
import csv
from pathlib import Path

# Define the directory of the current file
kfitdata_dir = Path('__file__').parent

# Output file path
output_file = 'kfitdata_combined.csv'

# Get all CSV files in the directory
csv_files = list(kfitdata_dir.glob('*.csv'))

if not csv_files:
    print("No CSV files found in kfitdata directory")
else:
    # Open output file for writing
    with open(output_file, 'w', newline='') as outfile:
        writer = None
        
        for csv_file in csv_files:
            # Get the filename stem (without extension)
            filename_stem = csv_file.stem
            
            with open(csv_file, 'r', newline='') as infile:
                reader = csv.DictReader(infile)
                
                # Initialize writer with headers on first file
                if writer is None:
                    fieldnames = ['Filename'] + reader.fieldnames
                    writer = csv.DictWriter(outfile, fieldnames=fieldnames)
                    writer.writeheader()
                
                # Write rows with filename stem added
                for row in reader:
                    row['Filename'] = filename_stem
                    writer.writerow(row)
    
    print(f"Combined {len(csv_files)} CSV files into {output_file}")