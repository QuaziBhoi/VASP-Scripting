#!/bin/bash

#Copy input files from static caluculation folder

cd ../dos
cp WAVECAR INCAR POTCAR KPOINTS ../ev/
cd ../ev
sed -i '/LWAVE/c\LWAVE = .FALSE.' INCAR

iteration=0

#For each lattice parameter in the specified length, generate a POSCAR file and run static calclation

for (( j=75; j<=125; j+=5 )); 
do
((iteration++))

echo "Iteration $iteration: Percentage = $j"
j_scaled=$(bc <<< "scale=16; $j / 100")
a=$(bc <<< "scale=16; 5.8012416593023870 * $j_scaled")
b=$(bc <<< "scale=16; -0.0142860073088748 * $j_scaled")
c=$(bc <<< "scale=16; -0.0361151569346229 * $j_scaled")
d=$(bc <<< "scale=16; 0.1356886841290714 * $j_scaled")
e=$(bc <<< "scale=16; 5.9852759620003972 * $j_scaled")
f=$(bc <<< "scale=16; 0.0166013975976411 * $j_scaled")
g=$(bc <<< "scale=16; 1.4635300798795077 * $j_scaled")
h=$(bc <<< "scale=16; 1.8609509820921870 * $j_scaled")
i=$(bc <<< "scale=16; 7.9015353820367862 * $j_scaled")
echo "a= $a" 
echo "b= $b"
echo "c= $c" 
echo "d= $d"
echo "e= $e" 
echo "f= $f"
echo "g= $g" 
echo "h= $h"
echo "i= $i"

cat >POSCAR <<! 
Bi2 Cr O6                               
   1.00000000000000     
     $a   $b    $c
     $d   $e    $f
     $g   $h    $i
   Bi   Cr   O 
     4     2    12
Direct
  0.9126381673016054  0.7614159078712374  0.1707298395287933
  0.0873618326983948  0.2385840921287552  0.8292701754712044
  0.5079040255091218  0.2705679787835618  0.1791065770390266
  0.4920959744908776  0.7294320212164451  0.8208934379609713
  0.6977550088496260  0.2568149835878814  0.5760261216538560
  0.3022449911503743  0.7431850164121253  0.4239738783461436
  0.8654132511110856  0.0485704245384291  0.6524580066195609
  0.1345867488889146  0.9514295714615776  0.3475419933804390
  0.4255761372485914  0.1544137967536235  0.6368610412020940
  0.5744238327514058  0.8455862182463745  0.3631389587979055
  0.7561050250863369  0.3440773181512080  0.3658653387880558
  0.2438949749136631  0.6559226518487967  0.6341346612119437
  0.6843085659399266  0.5958206514650328  0.0422509597578108
  0.3156914340600731  0.4041793485349668  0.9577490552421948
  0.7976456193131174  0.1007546107888134  0.0470748533771835
  0.2023543806868822  0.8992453822111893  0.9529251436228057
  0.7432678179348019  0.4647834508569781  0.6675890923683921
  0.2567321820651907  0.5352165491430213  0.3324109076316070
!

# Run VASP with the new POSCAR file
nohup mpirun -np 80 vasp 

# Gather energy and volume data from the VASP outputs and Store the data in "EvsV" text file

V=`grep volume/ion OUTCAR|awk '{print $5}'`
E=`tail -n1 OSZICAR | awk '{ print $5}'`
echo $V $E >> EvsV

# Clean up VASP files
rm OUTCAR OSZICAR

# Restart
  done
