# Colocalization

This folder contains 3 macros. The macro "Colocalization_Analysis.ijm" measures colocalization score for each region of interest and all channels. The macro "Local_Pearson_Coefficient.ijm" computes a map of Pearson Correlation Coefficient in 3D local Gaussian windows. The macro "Colocalization_By_Cell.ijm" is similar to the first one but here the ROI are defined by segmenting cells in the image from a nuclei labelling.

## Colocalization Analyzing

Compute various colocalization measures for all channels of an hyperstack for segmented structures.

If there is no opened image, a test image is generated.

![image](https://user-images.githubusercontent.com/3415561/117972406-05f37280-b323-11eb-970b-4adc9b49abb6.png)


### Tutorial
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


### Colocalization metrics
For each pair of channel, we consider the intensity above a channel-wise defined threshold for the preprocessed image. The colocalization metric can be performed on the original or on the pre-processed image enabling background substraction and denoising before colocalization. We denote x1 and x2 the list of intensity   preprocessed images in the union of objects in both channel and t1 and t2 the associated thresholds. The corresponding intensities in the original images are y1 and y2.

![image](https://user-images.githubusercontent.com/3415561/146014868-fc83c4e1-7e6c-4549-8250-5fe202e40ea1.png)


Here is the list of the 25 computed metrics:
- Pearson: [Pearson's correlation coefficient ](https://en.wikipedia.org/wiki/Pearson_correlation_coefficient) is a measure of linear correlation between two channels. It ranges from -1 (anticorrelation) to 1 (correlation).
- M1 : Manders' M1 score computed as sum in blue vs sum in blue and yellow. 
- M2 : Manders's M2 score computed as the sum in blue versus sum in blue and green
- MOC : Mander's overlap coefficient (Pearson without substracting the mean)
- D : Independance metric for the sets 
- N1 : number of segmented objects in channel 1
- N2 : number of segmented objects in channel 2
- N12 : number of overlapping objects 
- N12/N1 : ratio of overlapping objects over the objects in channel 1
- N12/N2 : ratio of overlapping objects over the objects in channel 2
- Area Object 1: Area of segmented objects in channel 1 in pixels
- Area Object 2: Area of segmented objects in channel 2 in pixels
- Area Intersection: Area of the intersection
- Intersection vs Area Object 1
- Intersection vs Area Object 2
- Mean Ch1 in Object 1
- Mean Ch2 in Object 2
- Mean Ch1 in Intersection
- Mean Ch2 in Intersection
- Mean Ch1 
- Stddev Ch1
- Max Ch1 : maximum on Ch1 (useful for detecting saturation)
- Mean Ch2
- Stddev Ch2
- Max Ch 2: maximum on Ch2 (useful for detecting saturation)


## Local Pearson Coefficient

This macro computes a map of Pearson Correlation Coefficient in 3D local Gaussian windows for two selected channels of a multi-channel hyperstack.

If no image is given a test image is generated. The image is showing two sin wave going out of phase along the horizontal axis.

![image](https://user-images.githubusercontent.com/3415561/122535614-a2dcb600-d01b-11eb-9767-e3e69d2015e6.png)

An the resulting map is:

![image](https://user-images.githubusercontent.com/3415561/122535701-b5ef8600-d01b-11eb-96e6-ad59cd51b7cd.png)

We can see the evolution of correlation coefficient going from 1 to -1 from left to right.


## Colocalization by Cell
The macro Colocalization_by_Cell.ijm is used to measure colocalization cell by cell in an image that would have cell nuclei labelled enabling an automated segmentation of the cells. It is similar to the macro 'Cell_Organization.ijm' as it can handle a single opened image, a file or a list of files for batch processing using a csv file as input.

Typing 'test' in the Filename parameter will generate a test image with 3 channels:

![image](https://user-images.githubusercontent.com/3415561/145996199-1493e912-8fda-440c-ba86-e2c69078b23d.png)

and run the macro with default parameters.

#### Format of the input csv file
The csv file contains a list of absolute or relative path to the files to process and optionally the first parameters of the macro in columns. The list of files can be generated using the macro 'File_Conversion/Parse_Folders.ijm'
- Filename  [compulsory] : absolute or relative path to the file
- Condition [optional] :a experimental group
- Channel Names [optional]: the comma-separated list of names for each channel in the stack, for example DAPI, Tfr, MyoVb
- ROI Channel [optional]: the index of the nuclei marker, for example 1
- Object Channel [optional]: the comma-separated list of channel indices used for colocalization, for example 2,3,4
- Mask Channel [optional]: index of an optional cell mask label or limiting the extent of cells using an additional channel.


## StartDist Overlap
This macro uses startdist 2D on every channel and planes and compute the overlap between the detected regions of interests

It requires a specific version of the startdist plugin which does not clear the ROI manager at each run ( https://github.com/jboulanger/stardist-imagej/releases/download/0.3.0-noclear/StarDist_-0.3.0.jar ). You can download the jar can copy it in place of the original file in th eplugin folder.

Open the macro and press run, define the name of the channel which will be used to report the results and the colors of the overlay that will help visualize the segmentation result. Both parameters are a comma separated list of characters:
![image](https://user-images.githubusercontent.com/3415561/167824587-d88d645a-9bce-455a-865d-951e884cc5c2.png)

The result is a table with a line per slice. If several images are processed the new lines are appended to the table.
![image](https://user-images.githubusercontent.com/3415561/167824727-d5c5d09a-cf7a-41a6-b9ae-13f9325ec3ce.png)

If there is no opened image, a test image is create:
![image](https://user-images.githubusercontent.com/3415561/167825894-5038eff1-06f7-4861-82a0-56596de69273.png)

Note that the macro does not amke any intensity measurement.


## Reference
 Manders E et al. (1993). "Measurement of colocalization of objects in dual-color confocal images." J Microsc Oxford 169:375–382.
 
 
 
