#!/bin/bash

### IMPORTANT ###
# Edit as required

### INSTRUCTIONS ###
# This script makes new directory "ftp" and runs finite temperature phonon calculation and data extraction
# Run this script using: bash finite_temp_phonon.sh

echo "-----------------------------------------------------------------------------------------------"

cd .. #Go to upper directory
mkdir alamode #Make new directory 'bs'
cd alamode #Enter directory 'bs'
mkdir 0_harmonic
cd 0_harmonic
cp  ../../relax/POSCAR ./ #Copy all files from 'dos' directory into current directory

total_atoms=$(awk 'NR==7 {sum=0; for(i=1;i<=NF;i++) sum+=$i; print sum}' POSCAR)
num_elements=$(awk 'NR==6 {print NF}' POSCAR)
elements=$(awk 'NR==6 {for(i=1;i<=NF;i++) printf "%s%s", $i, (i<NF ? " " : "")}' POSCAR)

# Reading values from POSCAR file
read -r a b c < <(sed -n '3p' POSCAR)
read -r d e f < <(sed -n '4p' POSCAR)
read -r g h i < <(sed -n '5p' POSCAR)


total_atoms=$(awk 'NR==7 {sum=0; for(i=1;i<=NF;i++) sum+=$i; print sum}' POSCAR) # Extract the total number of atoms from line 7

# Extract the elements and their counts from lines 6 and 7
elements=($(awk 'NR==6 {for(i=1;i<=NF;i++) print $i}' POSCAR))
counts=($(awk 'NR==7 {for(i=1;i<=NF;i++) print $i}' POSCAR))
current_element=1 # Initialize a variable to track the current element index
position_counter=0 # Initialize a variable to keep track of the current count of positions
last_line=$((9 + total_atoms - 1)) # Calculate the last line number for position data
output="Position Data with Element Index:\n" # Initialize a variable to store the output
while IFS= read -r line; do # Increment the position counter
    position_counter=$((position_counter + 1)) # Increment the position counter
    # Check if we need to move to the next element
    if (( position_counter > counts[current_element - 1] )); then
        current_element=$((current_element + 1))
        position_counter=1
    fi
    # Append the element index and position line to the output variable
    output+="$current_element $line\n"
done < <(awk "NR>=9 && NR<=$last_line" POSCAR)
output=${output%\\n} # Remove the trailing newline character
echo -e "$output" # Print the output variable

### Create INCAR ###
cat >alm.in1 <<!
  PREFIX = output
  MODE = suggest
  NAT = $total_atoms; NKD = $num_elements
  KD = $elements
/

&interaction
  NORDER = 1  # 1: harmonic, 2: cubic, ..
/

&cell
  1.88973 # factor in Bohr unit
  %a $b $c 
  $d $e $f 
  $g $h $i 
/

&cutoff 
  *-* None
/

&position
$output
/
!
