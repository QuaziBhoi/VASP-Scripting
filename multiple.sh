#!/bin/bash

### INSTRUCTIONS ###
# This script runs multiple bash files sequentially
# Run this script using:
# bash multiple.sh [script names sequentially]
# For example: 
# bash multiple.sh relax static 

# Check if any arguments are provided
if [ "$#" -eq 0 ]; then
    echo "No script names provided."
    exit 1
fi

# Loop through each argument (script name)
for script in "$@"; do
    # Check if the script file exists and is executable
    if [ ! -x "$script.sh" ]; then
        echo "Error: Script '$script' not found or not executable."
        exit 1
    fi
done

# Run each script sequentially
for script in "$@"; do
    echo "Running script: $script"
    source "$script.sh"
    if [ $? -ne 0 ]; then
        echo "Error: Script '$script' failed to execute properly."
        exit 1
    fi
done

echo "All scripts executed successfully."
