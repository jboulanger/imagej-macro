#@File (label="Filename",style="open") filename
#@String(label="Channels",value="all") channel
#@String(label="Slices",value="all") slice
#@String(label="Frames",value="all") frame
#@Boolean (label="Maximum Intensity Projection",value=false) mip
#@String (label="Mode",choices={"8-bit","16-bit","32-bit"}) mode
#@File (label="Output folder", style="directory") imageFolder
#@String(label="Tag", value="") tag
#@ String(label="Format",choices={"TIFF","PNG","JPEG"},style="list") format

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
print("Opening file " + filename);
open(filename);
id0 = getImageID();

if (!matches(channel, "all") || !matches(slice, "all") || !matches(frame, "all") ) {
	Stack.getDimensions(width, height, channels, slices, frames);
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

if (mip) {
	selectImage(id1);
	run("Z Project...", "projection=[Max Intensity]");
	id2 = getImageID();
	selectImage(id1); close();
	selectImage(id2);
}

run(mode);

name = File.getNameWithoutExtension(filename);
ext = getNewFileExtension(format);
oname = imageFolder + File.separator + name + tag + ext;
print("Saving image to file " + oname);
saveAs(format,oname);
close();

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