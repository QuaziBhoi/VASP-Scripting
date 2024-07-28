#!/bin/bash

### IMPORTANT ###
# Edit as required

### INSTRUCTIONS ###
# This script makes new directory "optical" and runs optical calculation
# Run this script using: bash optical.sh

echo "-----------------------------------------------------------------------------------------------"

cd ..
mkdir optical
cd optical
cp -r ../dos/* ./
cp INCAR INCAR.st

mm=$(grep "MAGMOM" INCAR.r | awk -F "=" '{print $2}')
nb=$((3*$(grep "number of bands    NBANDS" OUTCAR | awk -F "=" '{print $4}')))

mm_string="MAGMOM = $mm"
nb_string="NBANDS = $nb"

echo  -e  "101\nOP\n"  |  vaspkit > vaspkit.txt
sed -i '/NBANDS/d' INCAR #Delete pre-existing NBANDS tag in INCAR

sed -i "2 i $nb_string" INCAR #Insert NBANDS into INCAR
sed -i "2 i $mm_string" INCAR #Insert MAGMOM into INCAR
sed -i "2 i ISPIN = 2" INCAR #Insert ISPIN into INCAR

### Runing VASP ####
nohup mpirun -np 80 vasp

echo  -e  "711\n1\n"  |  vaspkit > vaspkit.txt

