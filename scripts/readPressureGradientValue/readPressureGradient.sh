#!/bin/bash

# Author: LB
# Date: 31.07.2019.
# Location: In the middle of nowhere between Hvar and Korcula.


inputFile=$1

# Check whether file is empty
if [ -z "$inputFile" ]; then
    echo "No log file selected!"
    exit
#Check whether file starts with "log." or ends with ".log"
elif [ ${inputFile: -4} = ".log" ] || [ ${inputFile:0:4} = "log." ]; then
    echo "Log file specified: $inputFile"
else
    echo "$inputFile is not a valid log file!"
    exit
fi

## Create file with only time step variables.

# Create column time in file time.dat
echo "Time" > time.dat

# Grep fetches all lines from inputFile starting with Time = ,
# awk prints the thrid word (Time = 300) which is a number.
grep -i '^Time = *' $inputFile | awk '{print $3}' >> time.dat

## Creates gradP files.

# Create first column
echo "Pressure gradient" > gradP.dat

# Grep fetches last line in each timestep which ALWAYS includes "ExecutionTime",
# prints this line and the three preceding lines (B3) which will include the
# last solution of the pressure gradient. Another grep command fetches only the
# part with pressure gradinet = some number and the last awk command only prints
# the last word aka the number, i.e. the pressure gradient value.


#grep -i -B 3 '^ExecutionTime = *' $inputFile | grep -i 'Pressure gradient = *' | awk '{print $NF}' >> gradP.dat

#EDIT: LB, 15/08/2019: Cut for fixed possitions with "," as delimiter. This is
#added so as to remove any value after the "Pressure Gradient =...," part of line.
#grep -i -B 3 '^ExecutionTime = *' $inputFile | grep -i 'Pressure gradient = *' | cut -d "," -f2 | awk '{print $NF}' >> gradP.dat

#EDIT: LB, 15/08/2019: Extyracts the part of the line with Pressure gradient and
# a number in a scientific notation. i.e. any number like 3.90890e-12
grep -i -B 3 '^ExecutionTime = *' $inputFile |  grep -oE '[Pp]ressure gradient = [0-9]{1,}.[0-9]{1,}e[+-][0-9]{1,}' | awk '{print $NF}' >> gradP.dat

#Combines time.dat and gradP.dat files into t_gradP.csv as two columns with ";"
# as the specified delimiter.
paste -d";" time.dat gradP.dat > t_gradP.csv

# Removes  garbage
rm time.dat gradP.dat
