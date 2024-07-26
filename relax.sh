#!/bin/bash

### IMPORTANT ###
# Edit as required - number of cores(80), ENCUT(500), MAGMOM(0), KPOINTS etc.

### INSTRUCTIONS ###
# This script makes new directory "relax" and runs structural relaxation
# Run this script using: bash relax.sh
# You can edit INCAR and KPOINTS as required inside the script by manually editing

echo "-----------------------------------------------------------------------------------------------"

# Check if POSCAR file exists
if [[ ! -f "POSCAR" ]]; then
  echo "POSCAR file not found!"
  exit 1
fi 

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

### Create INCAR ###
cat >INCAR <<!
Global Parameters
ISTART =  1            (Read existing wavefunction; if there)
ISPIN  =  2            (Non-Spin polarised DFT)
MAGMOM =  $MAGMOM      (Initial magnetic momentum)
LREAL  = .FALSE.       (Projection operators: automatic)
ENCUT  =  500          (Cut-off energy for plane wave basis set, in eV)
PREC   =  A            (Precision level)
LWAVE  = .TRUE.        (Write WAVECAR or not)
LCHARG = .TRUE.        (Write CHGCAR or not)
ADDGRID= .TRUE.        (Increase grid; helps GGA convergence)
KPAR   = $(($np/4))    (Divides k-grid into separate groups)
NCORE  = 8
 
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
    sleep 60  # Adjust the sleep duration as needed
done                        
