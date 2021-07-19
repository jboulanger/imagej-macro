#@File (label="Input file", style="open") filename
#@String (label="Channels", value="all", description="Use Duplicate formateg 1-3 or all") channel
#@String (label="Slices", value="all", description="Use Duplicate format eg 1-2 or all") slice
#@String (label="Frames", value="all", description="Use Duplicate format eg 1-10 or all") frame
#@String (label="Series", value="all", description="Extract series, use 'all' to extract all, 'last' to extract the last one") serie
#@Boolean (label="Maximum Intensity Projection", value=false, description="Perform a MIP on the image") mip
#@Boolean (label="Split channels", value=false, description="split channels") split_channels
#@String (label="Mode", choices={"Same","8-bit","16-bit","32-bit"}, description="Select the bit-depth of the image") mode
#@File (label="Output folder", style="directory") imageFolder
#@String (label="Format",choices={"TIFF","PNG","JPEG"},style="list") format
#@String(label="Tag", value="", description="Add a tag to the filename before the extension when saving the image") tag
#@boolean (label="Dummy run",value=true) dummy

/* Convert all series in a file to TIFs in a folder
 *  
 *  You can batch process files using the Batch function in the script editor
 *  
 *  Using the dummy mode enable to inspect image size as well and save the result in a cvs file
 *  
 * Jerome Boulanger for Marta 2021 - updated for Nathan
 */
 
run("Close All");
setBatchMode("hide");
name = File.getNameWithoutExtension(filename);
ext = getNewFileExtension(format);

run("Bio-Formats Macro Extensions");
print("Input file:" + filename);
Ext.setId(filename);
Ext.getSeriesCount(seriesCount);
print("File contains " + seriesCount + " series");

if (matches(serie,"all")) {
	s0 = 1;
	s1 = seriesCount
} if  (matches(serie,"last")) {
	s0 = seriesCount;
	s1 = seriesCount;
} else {
	str = split(serie,"-");
	s0 = parseInt(str[0]);
	s1 = parseInt(str[1]);
}


for (s = s0; s <= s1; s++) {
	
	str="open=["+filename+"] color_mode=Composite rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_"+s+"";	
	
	Ext.setSeries(s-1);
	Ext.getSeriesName(seriesName);
	Ext.getSizeX(sizeX);
	Ext.getSizeY(sizeY);
	Ext.getSizeZ(sizeZ);
	Ext.getSizeT(sizeT);
	if (dummy) {		
		print("bioformat args: " + str);	
		print("output file   : " + oname);		
		
	} else {
		
		print("Loading serie " + s + "/" + (s1-s0+1));			
		run("Bio-Formats Importer", str);

		id1 = processImage(channel, slice, frame, mip, mode);
		
		if (split_channels) {
			selectImage(id1);
			Stack.getDimensions(width, height, channels, slices, frames);
			for (c = 1; c <= channels; c++) {
				selectImage(
				run("Duplicate...", "duplicate channels="+c);
				idc = getImageID();
				oname = imageFolder + File.separator + name + "_serie_" + IJ.pad(s,4) + "_channel_" + IJ.pad(s,2) + tag + ext;	
				print("Saving serie "+ s +"/" + (s1-s0+1) + " to "+ oname);	
				saveAs(format,oname);	
				selectImage(idc);
				close();
			}
		} else {
			oname = imageFolder + File.separator + name + "_serie_" + IJ.pad(s,4) + tag + ext;
			print("Saving serie "+ s +"/" + (s1-s0+1) + " to "+ oname);	
			saveAs(format,oname);
		}
		selectImage(id1); close();		
	}
	setResult("Serie Index", s-1,s);
	setResult("Output  file", s-1, File.getName(oname));	
	setResult("Serie name", s-1, seriesName);
	setResult("Width", s-1, sizeX);
	setResult("Height", s-1, sizeY);
	setResult("Slices", s-1, sizeZ);
	setResult("Frames", s-1, sizeT);
	
}
Ext.close();
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

function processImage(channel, slice, frame, mip, mode) {
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
	
	if (!matches(mode,"Same")){
		run(mode);
	}
	
	return getImageID();
}