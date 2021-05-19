# GUV Intensity

These macros analyse the intensity of GUVs overtime. Two approaches are proposed. The first one is based on regions. The user define ROI to analyse and add them to the ROI manager. This allows to manually discard unwanted parts of the image where aggregates might be present and influcence the intensity measurements. The second approach is based on tracking the contours of the GUVs over time. The user initializes on the first frame an approximate contour of several selected GUVs. Aggregates are detected in a pre-processing step and discarded when measuring the average intensity on the contour.

When running these macros with no image opened, the macro generate a test image.



