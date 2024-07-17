#!/bin/bash
### starting interface ###
echo "--------------------------------------------------------------------------------"
echo "------------------------------------Welcome-------------------------------------"
echo "--------------------------------------------------------------------------------"
echo "------------For any queries please contact 'Quazi Shafayat Hossain'------------"
echo 
echo
### Setting VASP parameters ###
echo  -e  "111\n1\n"  |  vaspkit
echo  -e  "115\nBi\ns\nBi\np\nBi\nd\n"  |  vaspkit
cp PDOS_USER.dat Bi.dat
echo  -e  "115\nCr\ns\nCr\np\nCr\nd\n"  |  vaspkit
cp PDOS_USER.dat Cr.dat
echo  -e  "115\nO\ns\nO\np\n"  |  vaspkit
cp PDOS_USER.dat O.dat
rm PDOS_USER.dat

echo " Thanks for using the script. For any queries please contact 'Quazi Shafayat Hossain' "
