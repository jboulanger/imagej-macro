# File_Conversion

Two macros that help handling files in multi-serie container or a set of files.

The macros apply to only one file at thes using "Run". Use "Batch" to process multiple files.

## Export_Series.ijm
This macro exports the individual series from a container such as a LIF file acquired on a Leica microscope with the LAX software. It enable to apply a few common operations on channels, and stacks.

You can use the dummy run mode to simply export the name of the series to a table for example.

![image](https://user-images.githubusercontent.com/3415561/123614999-c9a9a200-d7fc-11eb-81af-565e10c09424.png)

Parameters:
- Input File: the file to process
- Channels : type all to keep all channels or the channels range (1-3) you want to keep as in the "Duplicate..." menu
- Slices : type all to keep all slices or the slices range you want to keep as in the "Duplicate..." menu
- Frames : type all to keep all frames or the frames range you want to keep as in the "Duplicate..." menu
- Series : type all to keep all or the range you want to keep
- Maximum Intensity Projection: Tick this box to compute the maximum intensty projection
- Mode : Select the mode (bit depth) of the images for saving. Select 'Same' to not preserve the original mode
- Output folder: Where files will be saved
- Format: the format of the files
- Tag: optinally add a tag to the file name.
- Dummy run: tick this box to test and not save files or just list the files in the serie.

The path of saved files are 'output folder' / 'image name'_serie_'indices'+'tag'.'extension'

## Convert_files.ijm

 This macro enable to convert files and apply transformation (MIP, extract channel, etc) on the fly.

![image](https://user-images.githubusercontent.com/3415561/123611958-088a2880-d7fa-11eb-8783-af43f79d9796.png)

Parameters:
- Input File: the file to process
- Channels : type all to keep all channels or the channels range (1-3) you want to keep as in the "Duplicate..." menu
- Slices : type all to keep all slices or the slices range you want to keep as in the "Duplicate..." menu
- Frames : type all to keep all frames or the frames range you want to keep as in the "Duplicate..." menu
- Series : type all to keep all or the range you want to keep
- Maximum Intensity Projection: Tick this box to compute the maximum intensty projection
- Mode : Select the mode (bit depth) of the images for saving. Select 'Same' to not preserve the original mode
- Output folder: Where files will be saved
- Format: the format of the files
- Tag: optinally add a tag to the file name.
- Dummy run: tick this box to test and not save files or just list the files in the serie.


The path of saved files are 'output folder' / 'image name'_serie_'indices'+'tag'.'extension'

## Parse_folder.ijm

This macro parses files in folder with 1 recursion depth and lists the files and top folders in a table. This table can be used as an input for other macros such as Convention/Cell_Organization.ijm and Cluster_Job/Launch_Job_Array.ijm


