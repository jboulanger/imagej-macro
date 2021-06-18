# Colocalization

This folder contains 2 macros. The macro "Colocalization_Analysis.ijm" measures colocalization score for each region of interest and all channels. The macro "Local_Pearson_Coefficient.ijm" computes a map of Pearson Correlation Coefficient in 3D local Gaussian windows.

## Colocalization Analyzing

Compute various colocalization measures for all channels of an hyperstack for segmented structures.

If there is no opened image, a test image is generated.

![image](https://user-images.githubusercontent.com/3415561/117972406-05f37280-b323-11eb-970b-4adc9b49abb6.png)


## Tutorial
1. Open an image
2. Add roi in the ROI manager (draw roi and press t)
3. Run the macro and set the parameters
   - Preprocessing: choose the type of preprocessing (Blobs (lapacian of Gaussian), Regions (median and rolling ball), None)
   - Scale : scale in pixel used in the preprocessing ste
   - False Alarm: control of the false alarm to set the threshold level
   - Area[px] : filter ROI by size
   - circularity: filter ROI by circularity
   - overlay color: list of colors for the overlay of each channel, for example: red green blue
   - channel tag: list of tags for each channel, for example: nuclei golgi mitochondria 
   - condition: add an experimental condition column to the table, for example: WT
   - Perform colocalization: allow to run the preprocessing without the colocalization step
   - Colocalize on original: colocalization statistics will be all computed on the original image otherwise the preprocessed image is used allowing to remove varying background for example. 
   - Check the results for each cell and inspect the overlays
4. Inspect the result table and the segmented overlay


## Local Pearson Coefficient

This macro computes a map of Pearson Correlation Coefficient in 3D local Gaussian windows for two selected channels of a multi-channel hyperstack.

If no image is given a test image is generated. The image is showing two sin wave going out of phase along the horizontal axis.

![image](https://user-images.githubusercontent.com/3415561/122535614-a2dcb600-d01b-11eb-9767-e3e69d2015e6.png)

An the resulting map is:

![image](https://user-images.githubusercontent.com/3415561/122535701-b5ef8600-d01b-11eb-96e6-ad59cd51b7cd.png)

We can see the evolution of correlation coefficient going from 1 to -1 from left to right.





