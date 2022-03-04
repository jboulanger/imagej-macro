# Align Channels

The macro 'Align_Channel.ijm' estimates the parameters of a 2D deformation model by detecting beads in a calibration image. Beads are detected and their coordinates are stored in a table. Then these coordinates are used to estimate the alignment of the channels using an iterative closest point algorithm. For this a set of basic linear algebra function was written to simplify matrix manipulations. Once estimated the transform can be stored into a file (Models.csv) and loaded to align an image. The alignement is done pixel by pixel in the ImageJ macro language.

To download the macro click [here](https://raw.githubusercontent.com/jboulanger/imagej-macro/main/Align_Channels/Align_Channels.ijm)

Limitations:
- The name of the model table has to be 'Models.csv'. Rename it if needed.
- The macro will work on small images (tested on 1024x1024) and is relatively slow.
- Some error may occur, try closing the results and models.csv tables.

## Tutorial
1. Open an hyperstack image of beads with multiple channels. If the stack is 3D, perform a maximum intensity projection.
If there are no image opened the macro generates a test image with 2 channels and beads.
![image](https://user-images.githubusercontent.com/3415561/117965271-777af300-b31a-11eb-84dd-4b5cb2a86bdd.png)

2. Run the macro and select the estimate mode to estimate the parameters

![image](https://user-images.githubusercontent.com/3415561/117965355-94afc180-b31a-11eb-84c6-d61f390d9592.png)

3. Save the table Models.csv 

![image](https://user-images.githubusercontent.com/3415561/117965418-a4c7a100-b31a-11eb-9531-cf4c28ec63be.png)

4. Open the image you want to align and the file Models.csv. The name of the table with parameters has to be 'Models.csv'

5. Run again the macro to align the image using the parameters previously estimated

![image](https://user-images.githubusercontent.com/3415561/117965474-b4df8080-b31a-11eb-9e87-021590724764.png)

6. Check the result

![image](https://user-images.githubusercontent.com/3415561/117965614-e5271f00-b31a-11eb-9ad2-09d67b6bdee1.png)


# Manually align channels
Instead of using beads to calibrate the transformation, this macro used ROI from the ROI manager to estimate the shift between images.

## Tutorial
1. Open an image and select matching ROIs in every channels and add then in the ROI manager (press the key 't'). Try spreading the ROI over the image as much as possible.
2. Run the macro in the estimate mode. Optionally save the table 'Models.csv'
3. Run the macro in the apply model with the Models.csv table opened.
4. Check the result


