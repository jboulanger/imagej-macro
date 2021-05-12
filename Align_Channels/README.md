# Align Channels

The macro 'Align_Channel.ijm' estimates the parameters of a 2D deformation model by detecting beads in a calibration image. Beads are detected and their coordinates are stored in a table. Then these coordinates are used to estimate the alignment of the channels using an iterative closest point algorithm. For this a set of basic linear algebra function was written to simplify matrix manipulations. Once estimated the transform can be stored into a file (Models.csv) and loaded to align an image. The alignement is done pixel by pixel in the ImageJ macro language.

If there are no image opened the macro generates a test image with 2 channels and beads.


## Tutorial
- Open an hyperstack image of beads with multiple channels. If the stack is 3D, perform a maximum intensity projection.
- Run the macro and select the estimate mode to estimate the parameters
- Save the table Models.csv 
- Open the image you want to align and the file Models.csv
- Run again the macro to align the image using the parameters previously estimated



