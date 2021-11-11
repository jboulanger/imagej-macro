#!/bin/tcsh
echo "Working directory: `pwd`"
echo "Hostname         : `hostname`"
echo "Starting time    : `date`"
# set the input file
set ifile = "$1"
# set the output file
set ofile = `echo $ifile | sed 's/\.tif/_processed\.tif/g'`
echo "Input file       : $ifile"
echo "Output file      : $ofile"
echo "Finishing time   : `date`"
