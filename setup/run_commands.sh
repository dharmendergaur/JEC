#!/bin/bash

# Define the MLTarget values
targets=("GenEt" "logGenEt" "GenEtByL1Et" "logGenEtByL1Et")
# Loop through the target values and run the command
for target in "${targets[@]}"; do
    nohup python3 calculate_L1JetSFs_usingBDT.py --MLTarget "$target" &> "out_${target}.txt" &
done

echo "All commands executed successfully."


