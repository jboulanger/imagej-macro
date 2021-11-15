#!/bin/tcsh
#SBATCH --mail-type=END,FAIL   # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --nodes=1              # Run all processes on a single node
#SBATCH --ntasks-per-node=1    #
#SBATCH --cpus-per-task=112    # Number of CPU cores per task
#SBATCH --time=01:00:00        # Time limit hrs:min:sec
echo "Working directory: `pwd`"
echo "Hostname         : `hostname`"
echo "Starting time    : `date`"
# set the input file
set ifile = "$1"
# set the output file
set ofile = `echo $ifile | sed 's/\.tif/_denoised\.tif/g'`
/lmb/home/jeromeb/bin/ndsafir -i $ifile -o $ofile
echo "Input file       : $ifile"
echo "Output file      : $ofile"
echo "Finishing time   : `date`"
