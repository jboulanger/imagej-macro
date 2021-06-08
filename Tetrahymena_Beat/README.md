# Tetrahymena beat

This set of tool estimate the beat cilia beating frequency of single tetrahymena imaged with transmitted light. Most of the issues are solved by tuning the acquisition so that a single specimen crawls within the field of view and with high enough frame rate so that the beating can be sampled correctly.

## Tutorial
- Run the macro in batch mode using the text editor in Fiji, this will enable you to select several files and process them at once.
- Save the table into a csv file
- Open the script process_kymographs.m in Matlab, adjust the path to the csv files and conditions names and run to estimate beating frequencies and produce graphs for all sequences. 
