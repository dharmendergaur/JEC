import os
import csv

def merge_csv_files(input_directory, output_file):
    with open(output_file, 'w', newline='') as out_csv:
        writer = csv.writer(out_csv)
        first_file = True
        files_merged = 0

        for file in os.listdir(input_directory):
            if file.startswith("L1Nano") and file.endswith(".csv"):  # Process only CSV files starting with L1Nano
                file_path = os.path.join(input_directory, file)
                with open(file_path, 'r', newline='') as in_csv:
                    reader = csv.reader(in_csv)
                    if first_file:
                        writer.writerow(next(reader))  # Write the header of the first file
                        first_file = False
                    else:
                        next(reader, None)  # Skip the header of subsequent files
                    for row in reader:
                        writer.writerow(row)
                files_merged += 1

    return files_merged

if __name__ == "__main__":
    input_directory = "/opt/ppd/scratch/dharmender_pjc96232/Root_files_JEC25rederivation/Muon0_jobs/out/CMSSW_14_0_7/src/"  # Change this to your directory path
    output_file = "output.csv"  # Output file name
    files_merged = merge_csv_files(input_directory, output_file)
    print(f"Merged CSV saved as {output_file}")
    print(f"Total files merged: {files_merged}")