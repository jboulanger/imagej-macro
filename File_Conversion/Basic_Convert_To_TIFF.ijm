#@File (label="Filename",style="open", description="Input filename, use Batch in the script editor to process multiple files") filename
#@File (label="Output folder", style="directory", description="Output folder where image will be saved") imageFolder
open(filename);
name = File.getNameWithoutExtension(filename);
oname = imageFolder + File.separator + name + ".tif";
saveAs("TIFF",oname);
close();