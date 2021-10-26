# Analysis of cell organization

This macro analyse the spread of fluorescent marker in cells by measuring features such the spread (2nd moment), entropy, area, distance to the nuclei.

![image](https://user-images.githubusercontent.com/3415561/136790160-5dd4153c-73cc-461b-bbd7-c9ea5210754f.png)


## Parameters
- Input file:
  - if an image is already opened, then this is ignored and the image is processed
  - single image file (tif,lsm,nd2,etc)
  - csv file with list of files to process with columns names Filename,Condition,Channel Names,ROI Channel,Object Channel, Mask Channel. Any other column will be reported in the tables as well.
  - if the name is ending in 'test', then a test image is processed
- Condition: condition or experiemental group. If the input file is a csv, the condition can be replaced by the one in the csv file.
- Channel names: comma separated list of name for each channel.  If the input file is a csv, the condition can be replaced by the one in the csv file.
- ROI channel: channel index  used to define ROI seeds (typically nuclei)
- Object channel: channel index used to define object of interest whose spread and intensity will be reported
- Mask channel: channel index used to mask background (typically a cell mask/membrane labelling). Set to 'max' to use a maximim intensity projection of all channels.
- Downsampling: images can be optionally downsampled for analysis
- Segment with StarDist: nuclei will be detected using startdist to define the ROI, otherwise a smooth/threshold/watershed approach is used
- ROI specificity: threshold used to detect ROI if not using stardist
- Object specificity: threshold used to detect object
- Remove ROI touching image borders: Tick this to remove object near the image border to avoid biases
- Make band: limit the analysis to a band around the ROI
- Mask background: retrict the analysis to foreground
- Distance min: minimum distance used for the histogram of distances
- Distance max: maximum distance used for the histogram of distances
- Additional measurements: comma separated list of measurements defined in the measurement tool. The values can be Area, Mean, StdDev, Mode, Min, Max, X, Y, XM, YM, Perim., BX, BY, Width, Height, Major, Minor, Angle, Circ., Feret, IntDen, Kurt, %Area, RawIntDen, Slice, FeretX,  FeretY,  FeretAngle, MinFeret, AR, Round, Solidity. The measurement will be applied to the object channel only. Mean intensity measurements for each channel are already included by default.

## Output
The macro creates 2 tables that can be save as comma separated (csv) files (can be opened in excel, Graphpad Prism, etc)
- Regions of Interest.csv : one line per ROI with Mean intensity, Spread, Entropy, etc per region of interest (cell).
- Objects.csv: one line per object with distance to ROI border (eg nuclei), mean intensity in each channel, and all 

The R script make_figures.R shows an example on how to load the file 'Regions of interest.csv' to produce boxplot per condition with pairwise t test using ggplot2.
