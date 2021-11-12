#!/bin/tcsh
#SBATCH --time=01:00:00
echo "Working directory: `pwd`"
echo "Hostname         : `hostname`"
echo "Starting time    : `date`"
echo "Array task id    : $SLURM_ARRAY_TASK_ID"
echo "Input file       : $1"
@ line = ( $SLURM_ARRAY_TASK_ID + 1 )
set row = `sed -n {$line}p $1`
echo "Row              : $row"
# edit here to set the relevant code
# "processing --filename $1 --index $SLURM_ARRAY_TASK_ID"
echo "Finishing time   : `date`"
