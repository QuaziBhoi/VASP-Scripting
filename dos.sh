#!/bin/bash

### IMPORTANT ###
# Edit as required

### INSTRUCTIONS ###
# This script makes new directory "dos" and runs static calculation, density of states data extraction
# Run this script using: bash dos.sh
# You can edit INCAR and KPOINTS as required inside the script by manually editing

echo "-----------------------------------------------------------------------------------------------"

cd .. #Go to upper directory
mkdir dos #Make new directory 'dos'
cd dos #Enter directory 'dos'
cp -r ../relax/* ./ #Copy all files from 'relax' directory into current directory

### Creating INCAR File ###
cp INCAR INCAR.r #Copy INCAR as INCAR.r for backup
echo  -e  "101\nST\n"  |  vaspkit > vaspkit.txt #Generate new INCAR using vaspkit
cp INCAR INCAR.st #Copy INCAR as INCAR.st for backup

line1=$(sed -n '/^$/=' INCAR.r | head -n 1) && echo $line1
line2=$(sed -n '/^$/=' INCAR.st | head -n 1) && echo $line2
sed -n '1,/^$/p' INCAR.r > INCAR
sed -n "$((line2 + 1)),\$p" INCAR.st >> INCAR
cp INCAR INCAR.st

#### Runing VASP ####
nohup mpirun -np 80 vasp

#### Generate DOS data ####
echo  -e  "111\n1\n"  |  vaspkit

elements=$(awk 'NR==6 {for (i=1; i<=NF; i++) print $i}' POSCAR | sort | uniq)
# Loop through each element
for element in $elements; do
    echo -e "115\n${element}\ns\n${element}\np\n${element}\nd\n" | vaspkit
    if [ -f PDOS_USER.dat ]; then
        cp PDOS_USER.dat "${element}.dat"
    else
        echo "PDOS_USER.dat not found for ${element}!"
    fi
done
rm -f PDOS_USER.dat

#Plot DOS
cp ~/bash/dosplot.py ./
python dosplot.py



