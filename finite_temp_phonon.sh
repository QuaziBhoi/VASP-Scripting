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
cp  ../../relax/* ./ #Copy all files from 'dos' directory into current directory

