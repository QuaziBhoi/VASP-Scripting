#!/bin/bash

### IMPORTANT ###
# Edit as required

### INSTRUCTIONS ###
# This script makes new directory "ftp" and runs finite temperature phonon calculation and data extraction
# Run this script using: bash finite_temp_phonon.sh

echo "-----------------------------------------------------------------------------------------------"

cd .. #Go to upper directory
mkdir ftp #Make new directory 'bs'
cd ftp #Enter directory 'bs'
cp -r ../dos/* ./ #Copy all files from 'dos' directory into current directory
