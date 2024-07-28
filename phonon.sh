#!/bin/bash

### IMPORTANT ###
# Edit as required

### INSTRUCTIONS ###
# This script makes new directory "phonon" and runs phonon calculation and data extraction
# Run this script using: bash phonon.sh
# # You can edit MAGMOM as required by manually editing i n this script

echo "-----------------------------------------------------------------------------------------------"

cd .. #Go to upper directory
mkdir phonon #Make new directory 'bs'
cd phonon #Enter directory 'bs'
cp ../relax/POSCAR ./ #Copy POSCAR from 'relax' directory into current directory

### Calling phonopy to construct cells ###
cp POSCAR POSCAR-unitcell
phonopy -d --dim="1 1 1" -c POSCAR-unitcell
cp SPOSCAR POSCAR

### Setting VASP parameters ###
echo  -e  "102\n2\n0.03\n"  |  vaspkit
read -p 'Number of processing core of your system: ' np

### Setting MAGMOM to 0 for all atoms ### You can edit MAGMOM as required 
IFS=' ' read -r -a counts <<< "$(sed -n '7p' POSCAR)"
MAGMOM=""
for count in "${counts[@]}"; do
  MAGMOM+="$count*0 "
done
MAGMOM=${MAGMOM% }  # Remove trailing space
echo $MAGMOM

cat >INCAR <<!
ISMEAR =  0            (Gaussian smearing)
SIGMA  =  0.05         (Smearing value in eV)
IBRION =  8            (determines the Hessian matrix using DFPT)
EDIFF  =  1E-08        (SCF energy convergence; in eV)
PREC   =  A            (Precision level)
ENCUT  =  500          (Cut-off energy for plane wave basis set, in eV)
IALGO  =  38           (Davidson block iteration scheme)
LREAL  = .FALSE.       (Projection operators: false)
LWAVE  = .FALSE.       (Write WAVECAR or not)
LCHARG = .FALSE.       (Write CHGCAR or not)
ADDGRID= .TRUE.        (Increase grid; helps GGA convergence)
#LEPSILON = T
#ISYM = 0 
KPAR  = 2
NSIM  = 8
ISPIN = 2
MAGMOM =  $MAGMOM 
!
#### Runing VASP ####
mpirun -np $np vasp

#### Runing VASP for NAC(non-analytical correction) ####
mkdir nac
cp POSCAR POTCAR KPOINTS nac/
cd nac

cat >INCAR <<!
ISMEAR =  0            (Gaussian smearing)
SIGMA  =  0.05         (Smearing value in eV)
IBRION =  -1           (determines the Hessian matrix using DFPT)
EDIFF  =  1E-08        (SCF energy convergence; in eV)
PREC   =  A            (Precision level)
ENCUT  =  500          (Cut-off energy for plane wave basis set, in eV)
IALGO  =  38           (Davidson block iteration scheme)
LREAL  = .FALSE.       (Projection operators: false)
LWAVE  = .FALSE.       (Write WAVECAR or not)
LCHARG = .FALSE.       (Write CHGCAR or not)
ADDGRID= .TRUE.        (Increase grid; helps GGA convergence)
LEPSILON = T
KPAR  = 2
#ISYM = 0 
NSIM  = 8
ISPIN = 2
MAGMOM =  $MAGMOM 
!
mpirun -np $np vasp
phonopy-vasp-born > BORN1
cp BORN1 BORN_Data
sed '/#/ c 14.399652' BORN1 > BORN
cp BORN ../
cd ../

phonopy --fc vasprun.xml

cat >python_auto.py <<!
import phonopy
from phonopy.phonon.band_structure import get_band_qpoints_and_path_connections

#path = [[[0, 0, 0], [0.5, 0, 0.5], [0.625, 0.25, 0.625]],
#        [[0.375, 0.375, 0.75], [0, 0, 0], [0.5, 0.5, 0.5], [0.5, 0.25, 0.75]]]
#labels = ["$\Gamma$", "X", "U", "K", "$\Gamma$", "L", "W"]
#qpoints, connections = get_band_qpoints_and_path_connections(path, npoints=51)
phonon = phonopy.load(unitcell_filename="POSCAR-unitcell", log_level=1,
                      supercell_matrix=[1, 1, 1],
                      primitive_matrix=[[1, 0, 0],
                                        [0, 1, 0],
                                        [0, 0, 1]])
#phonon.run_band_structure(qpoints, path_connections=connections)


phonon.auto_band_structure()

phonon.auto_band_structure(npoints=101, with_eigenvectors=False, with_group_velocities=False, plot=True, write_yaml=True).savefig( 'band_auto.png' , format = 'png', dpi = 400 )


phonon.auto_total_dos(mesh=100.0, is_time_reversal=True, is_mesh_symmetry=True, is_gamma_center=False, plot=True, write_dat=True, filename="total_dos.dat").savefig( 'dos_auto.png' , format = 'png', dpi = 400 )

phonon.plot_band_structure_and_dos().savefig( 'banddos.png' , format = 'png', dpi = 400 )

phonon.auto_projected_dos(mesh=100.0, is_time_reversal=True, is_gamma_center=False, plot=False, pdos_indices=False, legend=False, write_dat=True, filename="projected_dos.dat")

phonon.plot_projected_dos(pdos_indices=None, legend=None).savefig( 'pdos.png' , format = 'png', dpi = 400 )

phonon.run_thermal_properties(t_min=0, t_max=1000, t_step=10, temperatures=None, is_projection=False, band_indices=None, cutoff_frequency=None, pretend_real=False)

tprop_dict = phonon.get_thermal_properties_dict()

for t, free_energy, entropy, cv in zip(
        tprop_dict['temperatures'],
        tprop_dict['free_energy'],
        tprop_dict['entropy'],
        tprop_dict['heat_capacity']):
    print(("%12.3f " + "%15.7f" * 3) % (t, free_energy, entropy, cv), file=open("thermal.dat", "a"))
phonon.plot_thermal_properties().savefig ( 'thermal.png' , format = 'png', dpi = 400 )

!

python python_auto.py
phonopy-bandplot --gnuplot band.yaml > band.dat
cat >anime.conf <<!
DIM = 1 1 1
READ_FORCE_CONSTANTS = .TRUE.
ANIME_TYPE = XYZ
ANIME = 2 5 20 
# 2 is band index
# 5 is amplitude
# 20 in number of images
!
phonopy -p anime.conf
