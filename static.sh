#!/bin/bash
### starting interface ###
echo "--------------------------------------------------------------------------------"
echo

# read -p 'Number of processing core of your system: ' np

### Creating INCAR File ###
echo  -e  "101\nST\n"  |  vaspkit > vk.txt
cp INCAR INCAR.st

if grep -q "IVDW" INCAR.r; then
    sed -n '1,13p' INCAR.r > INCAR
else
    sed -n '1,12p' INCAR.r > INCAR
fi

sed -n '21,29p' INCAR.st >> INCAR

if grep -q "DFTU Calculation" INCAR.r; then
    if grep -q "IVDW" INCAR.r; then
    sed -n '28,36p' INCAR.r >> INCAR
    else
    sed -n '27,35p' INCAR.r >> INCAR
    fi
fi

#### Runing VASP ####
nohup mpirun -np 80 vasp

#### Generate DOS data ####
echo  -e  "111\n1\n"  |  vaspkit

echo  -e  "115\nBi\ns\nBi\np\nBi\nd\n"  |  vaspkit
cp PDOS_USER.dat PDOS_Bi
echo  -e  "115\nMo\ns\nMo\np\nMo\nd\n"  |  vaspkit
cp PDOS_USER.dat PDOS_Mo
echo  -e  "115\nO\ns\nO\np\n"  |  vaspkit
cp PDOS_USER.dat PDOS_O

#### Send Telegram Message ####
curl "https://api.telegram.org/bot6605743120:AAHsqkBpVtDt9vKtJRU7bcsSPLq7RLY3nx0/sendMessage?chat_id=6698838969&text={Finished Static Calculation}"

cp ~/bash/dosplot.py ./
python dosplot.py

curl -s -X POST "https://api.telegram.org/bot6605743120:AAHsqkBpVtDt9vKtJRU7bcsSPLq7RLY3nx0/sendDocument" \
-F "chat_id=6698838969" \
-F "document=@./DOS_Plot.png" \
-F "caption=DOS Plot"

