#!/bin/bash
### starting interface ###
echo "--------------------------------------------------------------------------------"

### Setting VASP parameters ###
cp INCAR INCAR.r
echo  -e  "101\nST\n"  |  vaspkit
cp INCAR INCAR.st
sed -i '1,21d' INCAR
cp INCAR INCAR.temp
sed -n '1,13p' INCAR.r > INCAR
sed -n '1,7p' INCAR.temp >> INCAR

# #### Runing VASP ####
nohup mpirun -np $np vasp

nb=$(grep "number of bands    NBANDS" OUTCAR | awk -F "=" '{print $4}')
mm=$(grep "MAGMOM" INCAR.r | awk -F "=" '{print $2}')
echo $nb
nbb=$((3*$nb))
echo $nbb

echo $mm
string="MAGMOM = $mm"
strig="NBANDS = $nbb"
echo $string
echo $strig

echo  -e  "101\nOP\n"  |  vaspkit > foo.txt
sed -i '3d' INCAR

sed -i "2 i $strig" INCAR
sed -i "2 i $string" INCAR
sed -i "2 i ISPIN = 2" INCAR

# #### Runing VASP ####
# nohup mpirun -np $np vasp

# echo  -e  "711\n1\n"  |  vaspkit

echo " Thanks for using the script. For any queries please contact 'Quazi Shafayat Hossain' "
