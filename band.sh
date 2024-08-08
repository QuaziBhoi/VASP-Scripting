#!/bin/bash

### IMPORTANT ###
# Edit as required

### INSTRUCTIONS ###
# This script makes new directory "bs" and runs band structure calculation and data extraction
# Run this script using: bash band.sh
# You can edit k-path as required by manually editing

echo "-----------------------------------------------------------------------------------------------"

cd .. #Go to upper directory
mkdir bs #Make new directory 'bs'
cd bs #Enter directory 'bs'
cp -r ../dos/* ./ #Copy all files from 'dos' directory into current directory

cp INCAR INCAR.st #Copy INCAR as INCAR.st for backup
sed -i '5 i ICHARG = 11' INCAR #Insert ICHARG tag into INCAR

#### Runing VASP ####
nohup mpirun -np 20 vasp_std

echo  -e  "303\n"  |  vaspkit > vaspkit.txt #Call vaspkit for auto generated k-path
cp KPOINTS KPOINTS.r #Copy KPOINTS as KPOINTS.r for backup
cp KPATH.in KPOINTS #Copy KPATH as KPOINTS for simulation

#### Runing VASP ####
nohup mpirun -np 20 vasp_std

echo  -e  "211\n1\n"  |  vaspkit > vaspkit.txt #Extract band structure data
cd ../_scripts
