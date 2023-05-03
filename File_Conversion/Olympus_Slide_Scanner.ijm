#@File    (label="Input file", description="Select a .vsi file", style="open") ipath
#@Integer (label="Resolution level", description="Select a resolution level to export the data (1-6) 1 is the biggest resolution", value=6) resolution_level
#@String  (label="Format", choices={"TIFF","PNG"}) format
#@String  (label="Channels",  description="List channels to export separated by a comma.Type 'all' for including all channels", value="all") channel_csv
#@String  (label="LUTs", description="List channels to export separated by a comma. Type 'Grays' for all channel in grayscale, look at the Image>Lookup Tables for other options.", values="auto") lut_csv
#@Boolean (label="Skip overview and labels", description="Skip overview and label images", value=true) skip_overview_label
#@String  (label="Summary", description="Kind of summary",choices={"Detailed","Short"}) summary
#@Boolean (label="Dummy run", description="Do not export files",value=false) dummy
#@File    (label="Output folder", style="directory") opath
#@String  (label="Tag", description="add a tag to the files",value="") tag

/*
 * Decode a multi resolution VSI file
 * 
 * PNG: use the first channel only
 * TIF: save all channels
 * 
 * Jerome for Elena and Ben
 */

default_channels = newArray(1,2,3,4,5,6,7,8,9,10);
channel_list = parseCSVInt(channel_csv, default_channels);

default_lut = newArray("auto");
lut_list = parseCSVString(lut_csv, default_lut);

info = parseSeriesNames(ipath);

if (matches(summary,"Detailed")) {
	displayInfo(info);
} else {
	displaySummary(info);
}

if (!dummy) {
	exportResolutionLevel(ipath, info, opath, resolution_level, skip_overview_label, format, channel_list, lut_list, tag);
}

run("Collect Garbage");


function exportResolutionLevel(ipath, info, opath, level, skip_overview_label, format, channel_list, lut_list, tag) {
	/*
	 * Export a selected resolution level
	 * 
	 * ipath : path to the original file
	 * info : information structure
	 * opath : output folder
	 * level : resolution level
	 * skip_overview : 1: will skip label, macri image and overview 
	 * format : file format for saving files (PNG or TIFF)
	 * channels_list : array with list of channels to keep
	 * lut_list : list of look up table to apply
	 * tag : tag to append to the files
	 */
	setBatchMode(true);	
	m = 10; // number of fiels in info
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
				selectChannels(channel_list, lut_list);
				if (channel_list.length > 1) {
					Stack.setDisplayMode("composite");				
					run("RGB Color");
				} else {
					run("8-bit");				
				}
				ofilename = basename + "_" + replace(name, ",", "_") +   tag + ".png";				
			} else {
				selectChannels(channel_list, lut_list);			
				if (channel_list.length > 1) {		
					Stack.setDisplayMode("composite");
				}
				ofilename = basename+ "_" + replace(name, ",", "_") + tag+ ".tif";
			}
			print("Saving " + ofilename);
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

function selectChannels(channel_list, lut_list) {
	/*	Select channels from a list and apply a LUT */	 
	id0 = getImageID();
	title = getTitle();
	Stack.getDimensions(width, height, channels, slices, frames);	
	if (channels > 1) {				
		run("Stack to Hyperstack...", "order=xyczt(default) channels="+channels+" slices="+slices+" frames=1 display=Color");	
		title = getTitle();
		
		run("Split Channels");
		wait(100);
		if (channel_list.length>1) {
			str = "";
			for (i = 0;	i < channel_list.length; i++) {
				ch = channel_list[i];
				if (ch <= channels) {										
					str += "c"+(ch)+"=[C"+(ch)+"-"+title+"] ";
				}
			}								
			run("Merge Channels...", str + " create keep");
		} else {
			run("Duplicate...", "title=["+title+"] duplicate ");
		}
		for (i = 0;	i < channel_list.length; i++) {
			Stack.setChannel(i+1);
			lut = lut_list[Math.min(i, lut_list.length-1)];				
			if (!matches(lut, "auto")) {
				run(lut_list[Math.min(i, lut_list.length-1)]);
			}
			run("Enhance Contrast...", "saturated=0.1");
		}		
		for (ch = 1; ch <= channels; ch++) {		
			selectWindow("C"+ch+"-"+title);
			close();
		}						
	}	
}

function parseCSVString(csv, default) {	
	if (matches(csv, "auto")) {
		values = default;
	} else {
		str = split(csv,",");
		values = newArray(str.length);
		for (i = 0 ; i < str.length; i++) {
			values[i] = String.trim(str[i]);
		}
	}
	return values;
}

function parseCSVInt(csv, default) {	
	if (matches(csv, "all")) {
		values = default;
	} else {
		str = split(csv,",");
		values = newArray(str.length);
		for (i = 0 ; i < str.length; i++) {
			values[i] = parseInt(str[i]);
		}
	}
	return values;
