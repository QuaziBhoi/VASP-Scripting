#!/bin/bash
### starting interface ###
echo "--------------------------------------------------------------------------------"
echo 

sed -i 's/ENCUT  =  500/ENCUT  =  620/g' INCAR
sed -i 's/IBRION =  2/IBRION =  6/g' INCAR
sed -i 's/KPAR/#KPAR/g' INCAR
sed -i 's/NCORE/#NCORE/g' INCAR

#### Runing VASP ####
nohup mpirun -np 20 vasp_std

echo  -e  "203\n"  |  vaspkit  >> Elastic_Data.txt
