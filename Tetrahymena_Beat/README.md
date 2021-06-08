# Tetrahymena beat

This set of tool estimate the beat cilia beating frequency of single tetrahymena imaged with transmitted light. Most of the issues are solved by tuning the acquisition so that a single specimen crawls within the field of view and with high enough frame rate so that the beating can be sampled correctly.

## Tutorial
- Open an 2D image sequence file in its original format (nd2) and run the Cilia_beat.ijm imageJ macro.
- Repeat the process for a few images
- Save the table into a csv file
- Open the script process_kymographs.m in Matlab, adjust the path, files and conditions and run to estimate beating frequencies and produce graphs for all sequences. 
