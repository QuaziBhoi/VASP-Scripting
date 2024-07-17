#!/bin/bash
### starting interface ###
echo "--------------------------------------------------------------------------------"
echo 

read -p 'Number of processing core of your system: ' np

mkdir elastic
cd elastic
cp -r ../relax/* ./
sed -i '/LWAVE/c\LWAVE = .FALSE.' INCAR
sed -i '/ENCUT/c\ENCUT =  620' INCAR
sed -i '/IBRION/c\IBRION =  6' INCAR
sed -i 's/NCORE/#NCORE/g' INCAR

#### Runing VASP ####
nohup mpirun -np 80 vasp

echo  -e  "203\n"  |  vaspkit  >> Elastic_Data.txt
