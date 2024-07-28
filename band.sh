#!/bin/bash

### IMPORTANT ###
# Edit as required

### INSTRUCTIONS ###
# This script makes new directory "bs" and runs band structure calculation and data extraction
# Run this script using: bash band.sh
# You can edit K=Path as required by manually editing

echo "-----------------------------------------------------------------------------------------------"

cd ..
mkdir bs
cd bs
cp -r ../dos/* ./

cp INCAR INCAR.st
sed -i '5 i ICHARG = 11' INCAR

#### Runing VASP ####
nohup mpirun -np 80 vasp

echo  -e  "303\n"  |  vaspkit > vaspkit.txt
cp KPOINTS KPOINTS.r
cp KPATH.in KPOINTS

nohup mpirun -np 80 vasp
echo  -e  "211\n1\n"  |  vaspkit > vaspkit.txt
