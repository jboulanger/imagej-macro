# GUV Intensity

These macros analyse the intensity of GUVs overtime. Two approaches are proposed. The first one is based on regions. The user define ROI to analyse and add them to the ROI manager. This allows to manually discard unwanted parts of the image where aggregates might be present and influcence the intensity measurements. The second approach is based on tracking the contours of the GUVs over time. The user initializes on the first frame an approximate contour of several selected GUVs. Aggregates are detected in a pre-processing step and discarded when measuring the average intensity on the contour.

When running these macros with no image opened, the macro generate a test image.


## Tutorial
### Using the region based approach
1. Load the image
2. Compute the maximum intensity projection
3. Select several ROI discarding the aggregates and add them to the ROI Manager.
4. Run the macro
5. Inspect the results and adjust parameters if necessary
6. Clear the ROI Manager and create new ROIs for estimating the background
7. Run the GUV_Background.ijm macro
8. Graph the correct signals using the macro GUV_Intensity_Graph.ijm 

### Using the contour based approach
1. Load the image
2. Create ROI on the 1st frame of the sequence close to the GUVs contours and add them to the ROI Manager
3. Run the macro
4. Inspect the results and adjust parameters if necessary
5. Optionnaly estimate and correct the background using the GUV_Background.ijm macro
6. Graph the correct signals using the macro GUV_Intensity_Graph.ijm 
