#@File (label="Input file", style="open") filename
#@String(label="channels to keep", description="Comma separated list of channel indices", value="1,2") keep_str
#@File (label="Output folder", style="directory") imageFolder
#@String (label="Format",choices={"TIFF","PNG","JPEG"},style="list") format
#@String(label="Tag", value="", description="Add a tag to the filename before the extension when saving the image") tag
/*
 *  Extract channels from an image stack
 *  
 *  Usage: press run and select the input file or run batch
 *  
 *  
 *  Jerome Boulanger 2021 for Marco Laub
 */

setBatchMode("hide");

// decode the coma separated channels and reorder them
keep = split(keep_str,",");
Array.sort(keep);

// open the file
print("Opening file " + filename);
str="open=["+filename+"] rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT";	
run("Bio-Formats Importer", str); 
run("Stack to Hyperstack...", "order=xyczt(default) channels="+nSlices+" slices=1 frames=1 display=Color");
name = File.getNameWithoutExtension(filename);
ext = getNewFileExtension(format);

// delete channels
Stack.getDimensions(width, height, channels, slices, frames);
for (c = channels; c > 0; c--) {
	deletechannel = true;
	for (i = 0; i < keep.length; i++) {
		if (keep[i]==c) {
			deletechannel = false;
		}
	}
	if (deletechannel) {
		Stack.setChannel(c);
		run("Delete Slice", "delete=channel");
	} else {
		run("Enhance Contrast", "saturated=0.35");
	}
}

run("Make Composite");

// save the image to a file and close
oname = imageFolder + File.separator + name + tag + ext;	
print("Saving image to file " + oname);
saveAs(format,oname);
close();
print("Done");

setBatchMode("exit and display");


function getNewFileExtension(format) {
	f = newArray("TIFF","PNG","JPEG");
	e = newArray(".tif",".png",".jpg");
	for (k = 0; k < f.length; k++) {
		if(matches(format, f[k])){
			k0 = k;
			return e[k];
		}
	}
	return ".tif";
}
