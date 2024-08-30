#!/bin/bash

### IMPORTANT ###
# Edit as required - number of cores(80), ENCUT(500), MAGMOM(0), KPOINTS etc.

### INSTRUCTIONS ###
# This script makes new directory "relax" and runs structural relaxation
# Run this script using: bash relax.sh
# You can edit INCAR and KPOINTS as required inside the script by manually editing

echo "-----------------------------------------------------------------------------------------------"

cd ..

# Check if POSCAR file exists
if [[ ! -f "POSCAR" ]]; then
  echo "POSCAR file not found!"
  exit 1
fi 

mkdir relax
cd relax
mv ../POSCAR ./
dos2unix POSCAR

### Setting VASP parameters ###
echo  -e  "102\n1\n0.02\n"  |  vaspkit > vaspkit.txt

### Setting MAGMOM to 0 for all atoms ### You can edit MAGMOM as required 
IFS=' ' read -r -a counts <<< "$(sed -n '7p' POSCAR)"
MAGMOM=""
for count in "${counts[@]}"; do
  MAGMOM+="$count*0 "
done
MAGMOM=${MAGMOM% }  # Remove trailing space
echo $MAGMOM

np=$(grep -c ^processor /proc/cpuinfo)
nc=$(echo "scale=0; sqrt($np)" | bc)
if [ "$nc" -ge 4 ]; then
    ncc=$nc
else
    ncc=4
fi

### Create INCAR ###
cat >INCAR <<!
Global Parameters
ISTART =  1            (Read existing wavefunction; if there)
ISPIN  =  2            (Non-Spin polarised DFT)
MAGMOM =  $MAGMOM      (Initial magnetic momentum)
LREAL  = .FALSE.       (Projection operators: automatic)
ENCUT  =  520          (Cut-off energy for plane wave basis set, in eV)
PREC   =  A            (Precision level)
LWAVE  = .TRUE.        (Write WAVECAR or not)
LCHARG = .TRUE.        (Write CHGCAR or not)
ADDGRID= .TRUE.        (Increase grid; helps GGA convergence)
KPAR   = $(($np / $ncc))             (Divides k-grid into separate groups)
NCORE  = $ncc
#IVDW   = 12
#LDAU    = .TRUE.        (Activate DFT+U)
#LDAUTYPE=  2            (Dudarev, only U-J matters)
#LDAUL   =  -1  2  -1    (Orbitals for each species)
#LDAUU   =  0  6  0      (U for each species)
#LDAUJ   =  0  0  0      (J for each species)
#LMAXMIX =  4            (Mixing cut-off, 4-d, 6-f)
#GGA = PS 
#GGA = CA

Electronic Relaxation
ISMEAR =  0            (Gaussian smearing; metals:1)
SIGMA  =  0.05         (Smearing value in eV; metals:0.2)
NELM   =  90           (Max electronic SCF steps)
NELMIN =  6            (Min electronic SCF steps)
EDIFF  =  1E-08        (SCF energy convergence; in eV)

Ionic Relaxation
NSW    =  100          (Max ionic steps)
IBRION =  2            (Algorithm: 0-MD; 1-Quasi-New; 2-CG)
ISIF   =  3            (Stress/relaxation: 2-Ions, 3-Shape/Ions/V, 4-Shape/Ions)
EDIFFG = -0.0001       (Ionic convergence; eV/AA)
!

echo "Please check the INCAR file for any errors or edits. Press any key except ENTER to open the file or wait 5 seconds to proceed automatically."
read -t 5 -n 1 user_input
if [ -n "$user_input" ]; then
  nano INCAR
else
  echo "No input received. Proceeding with the script."
fi

### Run VASP until reaching required accuracy ###
while true; do

    #### Runing VASP ####
    nohup mpirun -np $np vasp

    cp POSCAR POSCAR.old
    cp CONTCAR POSCAR
    cp INCAR INCAR.r

    # Check if the condition is met in the last two lines of nohup.out
    if tail -n 2 nohup.out | grep -q "reached required accuracy" || tail -n 2 nohup.out | grep -q "deleting file STOPCAR"; then
        break  # Exit the loop if the condition is met
    fi
    
    # Add a delay before the next iteration (optional)
    sleep 10  # Adjust the sleep duration as needed
done    
                    
sed -i '/IBRION/c\IBRION =  1' INCAR
nohup mpirun -np $np vasp
cp POSCAR POSCAR.old
cp CONTCAR POSCAR
cp INCAR INCAR.r

cd ../_scripts
