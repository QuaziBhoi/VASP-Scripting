#!/bin/bash

### IMPORTANT ###
# Edit as required

### INSTRUCTIONS ###
# This script makes new directory "elastic" and runs elastic calculation and data extraction
# Run this script using: bash elastic.sh

echo "-----------------------------------------------------------------------------------------------"

cd .. #Go to upper directory
mkdir elastic #Make new directory 'bs'
cd elastic #Enter directory 'bs'
cp -r ../relax/* ./ #Copy all files from 'dos' directory into current directory

sed -i '/LWAVE/c\LWAVE = .FALSE.' INCAR #Turn of wavefunction writing to save RAM
sed -i '/ENCUT/c\ENCUT =  620' INCAR #Change ENCUT to 620
sed -i '/IBRION/c\IBRION =  6' INCAR #Change IBRION to 6
sed -i 's/NCORE/#NCORE/g' INCAR #Turn off NCORE tag

#### Runing VASP ####
nohup mpirun -np 80 vasp

echo  -e  "203\n"  |  vaspkit  >> Elastic_Data.txt #Extract elastic data
