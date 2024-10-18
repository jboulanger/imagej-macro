# Quantification of calcium imaging signal

This macro aims at measuring the (F-F0)/F0 where F is the average fluorescence signal over time for predefine region of interest and F0 is the baseline.

- F0: the baseline is define as the a minimum filter of the median filtered signal
- SBR: signal to base ratio (F-F0)/F0
- zscore: normalized SRR by its variance SBR-mean_SBR/stddev_SBR
- peak: local maxima of the SBR or zscore above the threshold
- frequency: number of peaks over the duration of the sequence in mHz
- mean amplitude: mean ampltiude of the peaks

## Installation
Download the macro.

## Usage
1. Open images and draw ROIs and add them to the ROI manager
2. Save the ROI with the same name than the image but with the extensin _ROIset.zip or _ROIset.roi in the same folder
3. Close images, open the macro and press run to process a single file or batch to process several file.
4. Set the parameters

At the end of the process, you can double click on the file names to open the results files. The tables that are saved are:

- a table with the signals (-signal.csv). You can open this file and create graph in Fiji, by a right click on the table > Plot... Peak has a value corresponding to the amplitude of the peak.
  
   |Frame | roi-1-signal | roi-1-sbr | roi-1-zscore| roi-1-peak | roi-2-signal | roi-2-sbr | roi-2-zscore| roi-2-peak |  ... |
   |----|------|------|------|------|------|------|------|------|------|
   |0	|130.817901235|	0.026743538	|0.946918043|	0.000000000|	122.195121951|	0.020522483|	1.955992611|	0.000000000|...|
   |1	|130.817901235|	0.026743538	|0.946918043|	0.000000000|	122.195121951|	0.020522483|	1.955992611|	0.000000000|...|
  
  
- a table with the list of event (-events.csv)
  
  |roi | event | frame | time | amplitude |
  |----|------|------|------|------|
  |roi-1|0 |10|0.99|0.03|
  |roi-1|1 |10|0.99|0.03|
  |...|...|...|...|...|  
  |roi-2|2 |10|0.99|0.03|
  
- a table with statistics (-stats.csv)
  |roi | number | frequency | mean amplitude |
  |--|--|--|--|
  |roi-1|5|16|0.08|
  |roi-2|2|3|0.1|
  
- a figure with all signal (-figure.jpg)
