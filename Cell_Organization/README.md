# Analysis of cell organization

This macro analyse the spread of fluorescent marker in cells by measuring features such the 2nd moment, entropy, area, distance to the nuclei.




## Parameters
- Input file:
  - if an image is already opened, then this is ignored and the image is processed
  - single image file (tif,lsm,nd2,etc)
  - csv file with list of files to process with columns names Filename,Condition,Channel Names,ROI Channel,Object Channel. Any other column will be reported in the tables as well.
  - if the name is ending in 'test', then a test image is processed
- Condition: condition or experiemental group. If the input file is a csv, the condition can be replaced by the one in the csv file.
- Channel names: comma separated list of name for each channel.  If the input file is a csv, the condition can be replaced by the one in the csv file.
- Additional measurements: comma separated list of measurements defined in the measurement tool. The values can be Area, Mean, StdDev, Mode, Min, Max, X, Y, XM, YM, Perim., BX, BY, Width, Height, Major, Minor, Angle, Circ., Feret, IntDen, Kurt, %Area, RawIntDen, Slice, FeretX,  FeretY,  FeretAngle, MinFeret, AR, Round, Solidity. The measurement will be applied to the object channel only. Mean intensity measurements for each channel are already included by default.
