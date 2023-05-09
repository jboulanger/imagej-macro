#@File (label="Filename",style="open", description="Input filename, use Batch in the script editor to process multiple files") filename
#@String(label="Channels", value="all", description="Use Duplicate formateg 1-3 or all") channel
#@String(label="Slices", value="all", description="Use Duplicate formateg 1-2 or all") slice
#@String(label="Frames", value="all", description="Use Duplicate formateg 1-10 or all") frame
#@Boolean (label="Maximum Intensity Projection", value=false, description="Perform a MIP on the image") mip
#@String (label="Mode", choices={"Same","8-bit","16-bit","32-bit"}, description="Select the bit-depth of the image") mode
#@String (label="Display mode", choices={"composite","color","grayscale"}, description="Select display mode") display_mode
#@Boolean (label="Split Channels", value=false, description="Split each image in a different channel") split_channels
#@String(label="Contrast mode", choices={"Min/Max","Saturated 5%","Saturated 10%","None"}) contrast_mode
#@File (label="Output folder", style="directory", description="Output folder where image will be saved") imageFolder
#@String(label="Tag", value="", description="Add a tag to the filename before the extension when saving the image") tag
#@String(label="Format",choices={"TIFF","PNG","JPEG"}, style="list", description="File format as in the saveAs() command") format

/*
 * Convert files
 * 
 * Convert file enabling to select channel, slices, frames, maximum intensity projection
 * and conversion to 8/16/32 bits.
 * 
 * Jerome Boulanger for Nathan 2021
 * 
 */

setBatchMode(true);
print("- Opening file \n" + filename);
name = File.getNameWithoutExtension(filename);
ext = getNewFileExtension(format);

run("Bio-Formats Importer", "open=["+filename+"] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");

id1 = processStack(channel, slice, frame, mip, mode, display_mode, contrast_mode);
		
if (split_channels) {
	selectImage(id1);
	Stack.getDimensions(width, height, channels, slices, frame);
	for (c = 1; c <= channels; c++) {
		selectImage(id1);
		run("Duplicate...", "duplicate channels="+c);
		idc = getImageID();		
		oname = imageFolder + File.separator + name + "_channel_" + IJ.pad(s,2) + tag + ext;
		print("Saving channel" + c +" to \n"+ oname);	
		saveAs(format,oname);	
		selectImage(idc);
		close();
	}
} else {
	oname = imageFolder + File.separator + name + tag + ext;
	print("-Saving image to file\n"+ oname);	
	saveAs(format,oname);
}
selectImage(id1); close();

setBatchMode(false);

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

function processStack(channel, slice, frame, mip, mode, display_mode, contrast_mode) {
	/*
	 * Process a stack (Select planes, MIP, channel, bit depth)
	 * 
	 */
	
	id0 = getImageID();
	
	if (nSlices==1) {
		return id0;
	}
	
	Stack.getDimensions(width, height, channels, slices, frames);
	if (!matches(channel, "all") || !matches(slice, "all") || !matches(frame, "all") ) {
		if (matches(channel, "all")) {
			channel= "1-"+channels;
		}
		if (matches(slice, "all")) {
			slice = "1-"+slices;
		}
		if (matches(frame, "all")) {
			frame = "1-"+frames;
		}
		args = "channels="+channel+" slices="+slice+" frames="+frame;
		print(args);
		run("Duplicate...", "duplicate " + args);
		id1 = getImageID();
		selectImage(id0);close();
		selectImage(id1);
	} else {
		id1 = id0;
	}
	
	if (mip && slices > 1) {
		selectImage(id1);
		run("Z Project...", "projection=[Max Intensity]");
		id2 = getImageID();
		selectImage(id1); close();
		selectImage(id2);
	}
	
	for (channel= 1; channel <= channels; channel++) {
		Stack.setChannel(channel);
		if (matches(contrast_mode, "Min/Max")) {
			resetMinAndMax();
		} else if (matches(contrast_mode, "Saturated 5%")) {
			run("Enhance Contrast", "saturated=0.05");
		} else if (matches(contrast_mode, "Saturated 10%")) {
			run("Enhance Contrast", "saturated=0.1");
		} 
	}
	
	Stack.setDisplayMode(display_mode);
	
	if (!matches(mode,"Same")){
		run(mode);
	}
	
	return getImageID();
}