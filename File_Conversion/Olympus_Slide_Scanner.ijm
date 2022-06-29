#@File (label="Input file", style="open") ipath
#@Integer (label="Resolution level", value=1) resolution_level
#@String (label="Format", choices={"TIFF","PNG"}) format
#@Integer (label="Channel (png only)", value=4) channel_index
#@Boolean(label="Skip overview and labels", value=true) skip_overview_label
#@File (label="Output folder", style="directory") opath


/*
 * Decode a multi resolution VSI file
 * 
 * PNG: use the first channel only
 * TIF: save all channels
 * 
 * Jerome for Elena and Ben
 */

info = parseSeriesNames(ipath);
//displayInfo(info);
displaySummary(info);
exportResolutionLevel(ipath, info, opath, resolution_level, skip_overview_label, format);


function exportResolutionLevel(ipath, info, opath, level, skip_overview_label, format) {
	/*
	 * Export a selected resolution level
	 * ipath : path to the original file
	 * info : information structure
	 * opath : output folder
	 * level : resolution level
	 * skip_overview : 1: will skip label, macri image and overview 
	 * format : file format for saving files (PNG or TIFF)
	 */
	setBatchMode(true);	
	m = 10;
	n = info.length / m;
	k = 0;
	basename = File.getNameWithoutExtension(ipath);
	for (i = 0; i < n; i++) {
		name = info[m*i+1];
		if (info[m*i+2]==level && !( skip_overview_label && (matches(name,"label") || matches(name,"overview")))) {							
			run("Bio-Formats Importer", "open=["+ipath+"] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_"+(i+1));
			if (matches(format, "PNG")) {
				Stack.getDimensions(width, height, channels, slices, frames);
				if (slices > 1) {
					id0 = getImageID();
					run("Z Project...", "projection=[Max Intensity]");
					id = getImageID();
					selectImage(id0);close();
					selectImage(id);
				}
				Stack.setDisplayMode("grayscale");
				Stack.setChannel(channel_index);
				run("Enhance Contrast...", "saturated=0.1");
				run("8-bit");				
				ofilename = basename + "_" + IJ.pad(k+1, 3) + ".png";				
			} else {
				ofilename = basename+ "_" + IJ.pad(k+1, 3)+ ".tif";
			}			
			saveAs(format, opath + File.separator + ofilename);
			close();
			k++;		
		}
	}	
	setBatchMode(false);
	updateResults();
}

function getAllSeriesNames() {
	// Store all series names in an array
	Ext.getSeriesCount(n);
	names = newArray(n);
	for (i = 0; i < n; i++) {
		Ext.setSeries(i);
		Ext.getSeriesName(str);
		names[i] = str;
	}
	return names;
}

function displayInfo(info) {
	/*
	 * Display information for each serie (resolution level) in the vsi file
	 * info: array as 10 fields x Nseries
	 */
	run("Clear Results");
	m = 10;
	n = info.length / m;
	for (i = 0; i < n; i++) {
		setResult("Name",i,          info[m*i]);
		setResult("Type",i,          info[m*i+1]);
		setResult("Resolution",i,    info[m*i+2]);
		setResult("Magnification",i, info[m*i+3]);
		setResult("FOV",i,           info[m*i+4]);
		setResult("X",i,             info[m*i+5]);
		setResult("Y",i,             info[m*i+6]);
		setResult("Z",i,             info[m*i+7]);
		setResult("C",i,             info[m*i+8]);
		setResult("Index",i,         info[m*i+9]);
	}
	updateResults();
}

function displaySummary(info) {	
	/*
	 * Display information for each image in the vsi file
	 * info: array as 10 fields x Nseries
	 */
	run("Clear Results");
	m = 10;
	n = info.length / m;
	k = 0;
	for (i = 0; i < n; i++) {
		if (info[m*i+2]==1) {
			setResult("Name",k,          info[m*i]);
			setResult("Type",k,          info[m*i+1]);
			setResult("Resolution",k,    info[m*i+2]);
			setResult("Magnification",k, info[m*i+3]);
			setResult("FOV",k,           info[m*i+4]);
			setResult("X",k,             info[m*i+5]);
			setResult("Y",k,             info[m*i+6]);
			setResult("Z",k,             info[m*i+7]);
			setResult("C",k,             info[m*i+8]);
			setResult("Index",k,         info[m*i+9]);
			k++;
		}
	}
	updateResults();
}

function parseSeriesNames(path) {
	/*
	 * Open the image abd parse each
	 */
	run("Bio-Formats Macro Extensions");
	Ext.setId(path);
	// parse file names and return an info array with name,type, resolution level
	m = 10;
	Ext.getSeriesCount(n);	
	info = newArray(m*n);
	current_type = "";
	current_resolution = 1;	
	current_fov = "";
	current_mag = "";
	current_index = 0;	
	for (i = 0; i < n; i++) {
		Ext.setSeries(i);
		Ext.getSeriesName(name);
		Ext.getSizeX(sizeX);
		Ext.getSizeY(sizeY);
		Ext.getSizeZ(sizeZ);
		Ext.getSizeC(sizeC);
		Ext.getSizeT(sizeT);
		if (matches(name,"label")) {
			current_fov = "label";
			current_type = "label";
			current_resolution = 1;
		} else if (matches(name,"overview")) {
			current_fov = "overview";
			current_type = "overview";
			current_resolution = 1;
		} else if (matches(name,"macro image")) {
			current_fov = "macro image";
			current_type = "macro image";
			current_resolution = 1;
		} else if (matches(name, ".*vsi #.*")) {
			current_resolution++;
		} else if (matches(name,".*x_.*")) {
			parts = split(name,"_");			
			parts2 = split(parts[2],"#");
			current_mag = parts[0];
			current_fov = parseInt(parts2[0]);
			current_resolution = 1;
			current_type = name;
			current_index++;
		}
		info[m*i] = name;
		info[m*i+1] = current_type;
		info[m*i+2] = current_resolution;
		info[m*i+3] = current_mag;
		info[m*i+4] = current_fov;
		info[m*i+5] = sizeX;
		info[m*i+6] = sizeY;
		info[m*i+7] = sizeZ;
		info[m*i+8] = sizeC;	
		if (matches(current_type,"label") || matches(current_type,"overview") || matches(current_type,"macro image"))  {
			info[m*i+9] = "0";
		} else {
			info[m*i+9] = current_index;
		}
	}
	Ext.close();
	return info;
}

function openSerie(index) {	
	//TODO : maintain pixel size and colors..
	
	Ext.setSeries(index);
	Ext.getSeriesName(name);
	Ext.getSizeX(sizeX);
	Ext.getSizeY(sizeY);
	Ext.getSizeZ(sizeZ);
	Ext.getSizeC(sizeC);
	Ext.getSizeT(sizeT);
	//Ext.getPixelType();
	Ext.getImageCount(n);
	Ext.getDimensionOrder(dimOrder);
	newImage(name, "16-bit composite-mode", sizeX, sizeY, sizeC, sizeZ, sizeT);
	for (i = 0; i < n; i++) {
		Ext.getZCTCoords(i, z, c, t)
		Ext.openImage("plane", i);
		run("Copy");
		close();
		selectImage(id);
		Stack.setPosition(c, z, t);
		run("Paste");
	}
	id = getImageID();	
	Stack.setDisplayMode("composite");	
	id = getImageID();
	return id;
}

