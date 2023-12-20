# Parasite coat dynamic quantification

This macro quantifies the time evolution of the surface coating by a fluorescently tagged protein in 3D.
The first channel, labels the overall structure (inside and outside) of the parasite and the second channel 
on the 3D time lapse experiment labels the signal of interest on the surface.

![image](https://github.com/jboulanger/imagej-macro/assets/3415561/5533b6d2-4fee-4361-83d1-c999485a694d)

- Install MorpholibJ by enabling the update site Help>upadte>IJPB-plugins
- Open the image and crop the region you want to analyze with the object of interest centered on the first frame.
- Run the macro

An "Error function" model is fitted to the mean intensity to estimate the speed and time of the coating of the parasite.

![image](https://github.com/jboulanger/imagej-macro/assets/3415561/d5513fd7-c0ff-47f7-9063-a97775fac1cb)

The measurement are then exported to a table:

|Time mid point [min]|Time constant [min]|Intensity amplitude [au]|Intensity offset [au]|
|--------------------|-------------------|------------------------|---------------------|
|99.316548440|1.275861198|16.113436068|124.280343286|
