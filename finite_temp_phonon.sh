#!/bin/bash

### IMPORTANT ###
# Edit as required

### INSTRUCTIONS ###
# This script makes new directory "ftp" and runs finite temperature phonon calculation and data extraction
# Run this script using: bash finite_temp_phonon.sh

echo "-----------------------------------------------------------------------------------------------"

cd .. #Go to upper directory
mkdir alamode #Make new directory 'bs'
cd alamode #Enter directory 'bs'
mkdir 0_harmonic
cd 0_harmonic
cp  ../../relax/POSCAR ./ #Copy all files from 'dos' directory into current directory

total_atoms=$(awk 'NR==7 {sum=0; for(i=1;i<=NF;i++) sum+=$i; print sum}' POSCAR)
num_elements=$(awk 'NR==6 {print NF}' POSCAR)
elements=$(awk 'NR==6 {for(i=1;i<=NF;i++) printf "%s%s", $i, (i<NF ? " " : "")}' POSCAR)

# Reading values from POSCAR file
read -r a b c < <(sed -n '3p' POSCAR)
read -r d e f < <(sed -n '4p' POSCAR)
read -r g h i < <(sed -n '5p' POSCAR)

total_atoms=$(awk 'NR==7 {sum=0; for(i=1;i<=NF;i++) sum+=$i; print sum}' POSCAR) # Extract the total number of atoms from line 7

# Extract the elements and their counts from lines 6 and 7
elements=($(awk 'NR==6 {for(i=1;i<=NF;i++) print $i}' POSCAR))
counts=($(awk 'NR==7 {for(i=1;i<=NF;i++) print $i}' POSCAR))
current_element=1 # Initialize a variable to track the current element index
position_counter=0 # Initialize a variable to keep track of the current count of positions
last_line=$((9 + total_atoms - 1)) # Calculate the last line number for position data
while IFS= read -r line; do # Increment the position counter
    position_counter=$((position_counter + 1)) # Increment the position counter
    # Check if we need to move to the next element
    if (( position_counter > counts[current_element - 1] )); then
        current_element=$((current_element + 1))
        position_counter=1
    fi
    # Append the element index and position line to the output variable
    output+="$current_element $line\n"
done < <(awk "NR>=9 && NR<=$last_line" POSCAR)
output=${output%\\n} # Remove the trailing newline character

cat >alm.in1 <<!
&general
  PREFIX = output
  MODE = suggest
  NAT = $total_atoms; NKD = $num_elements
  KD = $(awk 'NR==6 {for(i=1;i<=NF;i++) printf "%s%s", $i, (i<NF ? " " : "")}' POSCAR)
/

&interaction
  NORDER = 1  # 1: harmonic, 2: cubic, ..
/

&cell
  1.88973 # factor in Bohr unit
  $a $b $c 
  $d $e $f 
  $g $h $i 
/

&cutoff 
  *-* None
/

&position
!

{
    echo -e "$output"
} >>alm.in1

{
    echo /
} >>alm.in1

alm alm.in1 > alm.log1

python ~/alamode/tools/displace.py --VASP=POSCAR --mag=0.01 -pf output.pattern_HARMONIC

num_files=$(ls disp*.POSCAR 2>/dev/null | wc -l) 
for ((i=1; i<=$num_files; i++))
do
    num=$(printf "%02d" $i)
    mkdir ${num}
    cd ${num}
    cp ../disp${num}.POSCAR POSCAR
    cp ../INCAR ./
    cp ../POTCAR ./
    cp ../KPOINTS ./
    nohup mpirun -np 80 vasp
    cp vasprun.xml vasprun${num}.xml 
    cp ./vasprun${num}.xml ../
    cd ../
done

python ~/alamode/tools/extract.py --VASP=POSCAR vasprun*.xml > DFSET_harmonic
cp alm.in1 alm.in2
sed -i 's/^  MODE = .*/  MODE = optimize/' alm.in2

sed -i '/&interaction/i\
&optimize\
  DFSET = DFSET_harmonic\
/' alm.in2

alm alm.in2 > alm.log2
grep "Fitting error" alm.log2

cat >phband.in <<!
&general
  PREFIX = output
  MODE = phonons
  FCSXML = output.xml
  NKD = $num_elements
  KD = $(awk 'NR==6 {for(i=1;i<=NF;i++) printf "%s%s", $i, (i<NF ? " " : "")}' POSCAR)
/

&cell
  1.88973 # factor in Bohr unit
  $a $b $c 
  $d $e $f 
  $g $h $i 
/

&kpoint
  1  # KPMODE = 1: line mode
!

#!/bin/bash

# Extract NPOINTS
npoints=$(grep "NPOINTS" KPATH.phonopy | awk '{print $3}')

# Extract BAND coordinates
band_coordinates=$(grep "BAND =" KPATH.phonopy | sed 's/BAND = //; s/,/ /g')

# Extract BAND labels
band_labels=$(grep "BAND_LABELS =" KPATH.phonopy | sed 's/BAND_LABELS = //; s/\$//g; s/\$//g')

# Convert the BAND coordinates into an array
read -a coords <<< "$band_coordinates"
# Convert the BAND labels into an array
read -a labels <<< "$band_labels"

# Initialize an empty array for the formatted K-path
kpath=()

# Loop through the coordinates and labels to format the K-path
for ((i=0; i<${#labels[@]}-1; i++)); do
    # Get the current and next label
    label1=${labels[i]}
    label2=${labels[i+1]}
    
    # Get the coordinates for the current and next K-point
    x1=${coords[i*3]}
    y1=${coords[i*3+1]}
    z1=${coords[i*3+2]}
    
    x2=${coords[(i+1)*3]}
    y2=${coords[(i+1)*3+1]}
    z2=${coords[(i+1)*3+2]}
    
    # Format the line and add it to the K-path array
    kpath+=("  $label1 $x1 $y1 $z1 $label2 $x2 $y2 $z2 $npoints")
done

# Append the formatted K-path to the file phband.in
{
    printf "%s\n" "${kpath[@]}"
} >> phband.in

{
    echo /
} >> phband.in

###Anharmonic Force

