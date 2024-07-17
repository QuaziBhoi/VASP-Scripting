#!/bin/bash

read -p 'Number of processing core of your system: ' np

cp -r ../dos/* ./
cp INCAR INCAR.st

sed -i '5 i ICHARG = 11' INCAR

nohup mpirun -np $np vasp

echo  -e  "303\n"  |  vaspkit > vaspkit.txt
cp KPOINTS KPOINTS.r
cp KPATH.in KPOINTS


nohup mpirun -np $np vasp
echo  -e  "211\n1\n"  |  vaspkit > vaspkit.txt
