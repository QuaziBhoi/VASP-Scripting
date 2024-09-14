#!/bin/bash

### IMPORTANT ###
# Edit as required - number of cores(80), ENCUT(500), MAGMOM(0), KPOINTS etc.

### INSTRUCTIONS ###
# This script makes new directory "relax" and runs structural relaxation
# Run this script using: bash relax.sh
# You can edit INCAR and KPOINTS as required inside the script by manually editing

echo "-----------------------------------------------------------------------------------------------"

cd ..

#Copy input files from static caluculation folder
mkdir ev
cd ev
cp ../dos/{POSCAR,WAVECAR,INCAR,POTCAR,KPOINTS} ./ 
cp POSCAR POSCAR.r
sed -i '/LWAVE/c\LWAVE = .FALSE.' INCAR

# Reading values from POSCAR file
read -r a b c < <(sed -n '3p' POSCAR)
read -r d e f < <(sed -n '4p' POSCAR)
read -r g h i < <(sed -n '5p' POSCAR)

iteration=0
#For each lattice parameter in the specified length
#Generate a POSCAR file and run static calclation
for (( j=75; j<=125; j+=5 )); # 75% to 125%
do
((iteration++))

echo "Iteration $iteration: Percentage = $j"
j_scaled=$(bc <<< "scale=16; $j / 100")
a_scaled=$(bc <<< "scale=16; $a * $j_scaled")
b_scaled=$(bc <<< "scale=16; $b * $j_scaled")
c_scaled=$(bc <<< "scale=16; $c * $j_scaled")
d_scaled=$(bc <<< "scale=16; $d * $j_scaled")
e_scaled=$(bc <<< "scale=16; $e * $j_scaled")
f_scaled=$(bc <<< "scale=16; $f * $j_scaled")
g_scaled=$(bc <<< "scale=16; $g * $j_scaled")
h_scaled=$(bc <<< "scale=16; $h * $j_scaled")
i_scaled=$(bc <<< "scale=16; $i * $j_scaled")

sed -i "3s/.*/  $a_scaled $b_scaled $c_scaled/" POSCAR
sed -i "4s/.*/  $d_scaled $e_scaled $f_scaled/" POSCAR
sed -i "5s/.*/  $g_scaled $h_scaled $i_scaled/" POSCAR

# Run VASP with the new POSCAR file
nohup mpirun -np $(grep -c ^processor /proc/cpuinfo) vasp

# Gather energy and volume data from the VASP outputs
# Store the data in "EvsV" text files
V=`grep 'volume of cell' OUTCAR | awk '{print $5}' | head -n 1`
E=`grep 'free  energy   TOTEN  =' OUTCAR | awk '{print $5}'`
echo $V $E >> EvsV

# Clean up VASP files
rm OUTCAR OSZICAR

# Restart
  done


python ../_scripts/evfit.py

