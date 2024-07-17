#!/bin/bash
### starting interface ###
echo "--------------------------------------------------------------------------------"
echo "------------------------------------Welcome-------------------------------------"
echo "--------------------------------------------------------------------------------"
echo "------------For any queries please contact 'Quazi Shafayat Hossain'-------------"
echo 
echo
while true; do
    read -p "
Please make sure there are only the 'POSCAR' and the 'relax.sh' files in this working directory.

If the above mentioned conditions are not satisfied then please enter 'no' and modify your working directory and run this bash file again. If it is yes then enter 'yes' to continue: " yn
echo
    case $yn in
        [Yy]* ) break;;
        [Nn]* ) exit;;
        * ) echo "Please answer yes or no.";;
    esac
done

read -p 'Number of processing core of your system: ' np
read -p 'Energy cutoff: ' ec


## Setting VASP parameters ###
echo  -e  "102\n1\n0.02\n"  |  vaspkit > vaspkit.txt

cat >INCAR <<!
Global Parameters
ISTART =  1            (Read existing wavefunction; if there)
ISPIN  =  2            (Non-Spin polarised DFT)
MAGMOM = 8*2 4*2 4*2 24*2
LREAL  = .FALSE.          (Projection operators: automatic)
ENCUT  =  $ec          (Cut-off energy for plane wave basis set, in eV)
PREC   =  A            (Precision level)
LWAVE  = .TRUE.        (Write WAVECAR or not)
LCHARG = .TRUE.        (Write CHGCAR or not)
ADDGRID= .TRUE.        (Increase grid; helps GGA convergence)
KPAR   = $(($np/4))            (Divides k-grid into separate groups)
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
