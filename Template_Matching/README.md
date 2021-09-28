# Template Matching

This folder contains an implementation of template matching as an ImageJ macro and and example of use for detecting large microsphere in 2D confocal images.

## Template_Matching.ijm
This macro identify similar objects than the one in the ROI manager in the current image using normalized cross-correlation. 

To test the macro, run it with no other image opened. It will use the sample image Dot Blot.jpg to and init the ROI Manager with a ROI located on the top left most dot.

![image](https://user-images.githubusercontent.com/3415561/135070992-1f0a28ae-d47f-464c-8dc7-5f92b3df936d.png)

The parameters of the macro are
- channel : the channel on which the correlation is computed
- downsampling : the downsampling factor used to speed up computation
- threshold : the threshold used to identified peaks. Perfect correlation is 1.
- average : If ticked, then an average of the matched template will be performed here. Templates are not registered before averaging.

## Bead_Intensity.ijm

This macro uses cross correlation or dark blob detection to identify large micro-sphere in the image.
Once localized, a contour based approach is used to localize the contour with subpixel accuracy and measure intensities in every channel.

To test the macro, run it with no opened image. A test image with be created. 
![image](https://user-images.githubusercontent.com/3415561/135072233-53b1716d-9903-4ebb-9dba-cc21f8b19a73.png)

The parameters are the following:
- Channel: the index of the channel used to localized the micro-spheres
- Channel names: the name of the channels used to report the intensity
- Mode : type of detection mode (Correlation, Blobs, skip), if ROI are already present or manually selected use skip
- Downsampling : downsampling factor used to speed up computation for correlation
- Treshold : the threshold (-1,1) value used for the correlation approach
- Diameter : the typical diameter in um of the objects. This is used to filter ROI after fitting
- Smoothness : smoothness of the curves, the higher the more round and less wiggely they'll look
- Bandwidth [px]: the with of the band that is used to measure intensity

