# FRAP measure

This macro helps measuring fluorescence recovery after photoactivation (FRAP) when the sample undergo some motion. Since the apprearance of region of interest will change dramatically after bleaching tracking can become unreliable, hence we use here a simple set of control points over time to interpolate the trajectory. The ijm file contains 3 macros with 3 hotkeys:
- [z] find the brightest plane in the z stack, 
- [g] refine the location of the brightest structure in x and y
- [q] measure intensties along the track

These can be combine with the ROI Manager's hotkey: [t] and [u] to define and refine a list of control points. To properly use this macro, it needs to be installed in order to activate the hotkeys. Several tracks can be defined in parallel and the macro will cluster the control points to identify the different tracks as long as these control point remains at a reasonable distance from each other. This allows to measure simulateanously the FRAP recovery, the background level by choosing a dark region and a control region to estimate bleaching. Finally, background and bleaching correction allow to normalize the FRAP intensity and further estimate the parameter of an exponential recovery, providing some insights on the recovery rates. Measured intensity can be save to a csv file to be later analyzed in depth.



## Tutorial
1. Save the macro and install it in ImageJ (plugin>macro>install)
2. Open you image with bioformat importer (to get the time stamps correct)
3. For each region of FRAP, Control (in cell) and Background:
   - Move the time slider to the 1st time point
   - Create a ROI (circle), the size is not too important
   - Optimize its location in xy by pressing the key 'g' and in z by pressing the key 'z'
   - Add the ROI to the ROI manager
   - Move to another time point (as long the movement is linear between time point it is ok), and repeat until you reach the end of the movie
4. Press the key 'q' to quantify the FRAP recovery
5. Assign FRAP, Control and Background labels to the ROIs and set the time of photobleaching (in frame).

# Morph_ROI
This macro helps measure intensity in ROIs along time. From a define set of ROI, the macro interpolates the shapes of ROI along time to populate the missing frames using the existing ROI as control points. Several group of ROI can be defined by giving them the same name.

## Tutorial
1. Download and save the [macro](https://raw.githubusercontent.com/jboulanger/imagej-macro/main/FRAP_Measure/Morph_ROI.ijm)
2. Open an image sequence
3. Draw shapes at various time points and add them to the ROI manager
4. Using the "rename" tool of the ROI manager, define groups of the ROIs (eg cell and foci)
5. Run the macro and list the name of the channels using comma as separator, list the name of ROI you want to use and their colors. Finally choose a measurement in the list.
6. Press run

Running the macro with no opened image will create a test image and should produce this montage:

![image](https://user-images.githubusercontent.com/3415561/184339037-7b90ced4-c9e9-4b6a-834c-ee6b32a30bc0.png)

along with those measutements:

![image](https://user-images.githubusercontent.com/3415561/184339208-e6e5b23a-fb6f-46c5-baac-4d8c9fa67611.png)

