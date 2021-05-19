# GUV Intensity

These macro analyse the intensity of GUVs overtime. Two approaches are proposed. The first one is based on regions. The user define regions to analyse and add them to the ROI manager. This allows to manually discard unwanted regions of the image where aggregate might be present and influcence the intensity measurements. The second approach is based on tracking the contours of the GUVs over time. The user initialize on the first frame using the circle selection for example the initial approximate contour of several GUVs. Here aggreate are detected in a pre-processing step and discarded when measuring the average intensity on the contour.

