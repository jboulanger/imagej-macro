#@File(label="Input file") filename
#@Float(label="Pixel size", value = 0.062) pixel_size
#@String(label="Unit", value="micron") unit
#@String(label="Pixel type", choices={"8-bit","16-bit","32-bit"}) pixel_type

/*
 * Set Pixel Size
 * 
 * Set the pixel size and type for the image and overwrite the image.
 * 
 * Warning this will update metadata.
 * 
 */
 
setBatchMode(true);
open(filename);
run("Set Scale...", "distance=1 known="+pixel_size+" unit="+unit);
run(pixel_type);
saveAs("TIFF", filename);
print(filename);
close();
setBatchMode(false);