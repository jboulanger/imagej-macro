#@File (label="Input file", style="open") filename
#@String (label="Channels", value="all", description="Use 'Make Subset' format eg 1-3 / 1,4,6 or all") channel
#@String (label="Slices", value="all", description="Use 'Make Subset' format eg 1-2 or all") slice
#@String (label="Frames", value="all", description="Use 'Make Subset' format eg 1-10 or all") frame
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
 * Jerome Boulanger for Marta 2021, updated for Nathan
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
	s1 = seriesCount;
} else if  (matches(serie,"last")) {
	s0 = seriesCount;
	s1 = seriesCount;
} else {
	str = split(serie,"-");
	s0 = parseInt(str[0]);
	s1 = parseInt(str[1]);
}

nR0 = nResults;
for (s = s0; s <= s1; s++) {
	
	str="open=["+filename+"] color_mode=Composite rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_"+s+"";	
	
	Ext.setSeries(s-1);
	Ext.getSeriesName(seriesName);
	Ext.getSizeX(sizeX);
	Ext.getSizeY(sizeY);
	Ext.getSizeZ(sizeZ);
	Ext.getSizeC(sizeC);
	Ext.getSizeT(sizeT);
		
	print("Loading serie " + s + "/" + (s1-s0+1));			
	if (!dummy) {
		run("Bio-Formats Importer", str); 
		id1 = processImage(channel, slice, frame, mip, mode);
	}
	
	if (split_channels) {
		for (c = 1; c <= sizeC; c++) {
			oname = imageFolder + File.separator + name + "_serie_" + IJ.pad(s,4) + "_channel_" + IJ.pad(c,2) + tag + ext;	
			if (!dummy) {
				selectImage(id1); 				
				run("Duplicate...", "duplicate channels="+c);
				idc = getImageID();								
				print("Saving serie "+ s +"/" + (s1-s0+1) + " to "+ oname);	
				saveAs(format,oname);	
				selectImage(idc);
				close();
			} else {
				print("Serie " + s + "Channel "+c+ " will be saved in " + oname);
			}			
		}
	} else {
		oname = imageFolder + File.separator + name + "_serie_" + IJ.pad(s,4) + tag + ext;
		if (!dummy) {
			print("Saving serie "+ s +"/" + (s1-s0+1) + " to "+ oname);	
			saveAs(format,oname);
		} else {
			print("Serie " + s + " will be saved in " + oname);
		}
	}
	
	if (!dummy) {selectImage(id1); close();	}

	setResult("Serie Index", nR0, s);
	setResult("Output  file", nR0, File.getName(oname));	
	setResult("Serie name", nR0, seriesName);
	setResult("Width", nR0, sizeX);
	setResult("Height", nR0, sizeY);
	setResult("Slices", nR0, sizeZ);
	setResult("Channels", nR0, sizeC);
	setResult("Frames", nR0, sizeT);
	nR0 = nR0 + 1;
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
		run("Make Subset...",args");		
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
