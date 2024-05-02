//@File    (label="Input file", description="Select a .vsi file", style="open") ipath
//@String  (label="Action", choices={"Show information","Export Images","Process Images"}) action
//@Integer (label="Resolution level", description="Select a resolution level to export the data (1-6) 1 is the biggest resolution", value=6) resolution_level
//@String  (label="Look up table", description="comma separated LUT (Blue,Green,Red,Grays, etc) or dyes (DAPI,FITC,Cy3,Cy5,Cy7) or mix") luts
//@String  (label="Subset", value="all",description="subset of image eg: channels=1-4 slices=1-5 or all") subset
//@String  (label="Projection", choices={"None","Average Intensity","Max Intensity","Min Intensity","Median"}) projection
//@File    (label="Output folder", style="directory") opath
//@String  (label="Tag", description="add a tag to the files",value="") tag
//@String  (label="Format", choices={"TIFF","PNG","JPEG"}) format
//@Integer (label="Tile", value=100) tile_size
//@File    (label="Macro",value="None",description="Name of the macro for processing images") macroname

/*
 * Decode, export, process a multi resolution VSI file
 *  
 * Jerome for Elena and Ben
 */

function parseCellSenseSeries(path) {
	/*
	 * Parse file and return an info array with 13 fields:
	 * 0: name,
	 * 1: type label/overview/image/macroimage
	 * 2: resolution level (1-6)
	 * 3: pixel size xy
	 * 4: pixel size z
	 * 5: objective mag
	 * 6: field of view
	 * 7: channel names
	 * 8: width
	 * 9: height
	 * 10: depth
	 * 11: number of channels
	 * 12: image index
	 */
	m = 13; // number of fields
	Ext.setId(path);
	Ext.getSeriesCount(n);
	info = newArray(m*n);
	current_type = "";
	current_resolution = 1;
	current_fov = "";
	current_mag = "";
	current_index = 0;
	current_full_sizeX = 1;
	current_full_sizeZ = 1;
	
	for (i = 0; i < n; i++) {
		Ext.setSeries(i);
		Ext.getSeriesName(name);
		Ext.getSizeX(sizeX);
		Ext.getSizeY(sizeY);
		Ext.getSizeZ(sizeZ);
		Ext.getSizeC(sizeC);
		Ext.getSizeT(sizeT);
		Ext.getPixelsPhysicalSizeX(current_pixel_size_x);		
		Ext.getPixelsPhysicalSizeZ(current_pixel_size_z);
		if (matches(name,"label")) {
			current_fov = "label";
			current_type = "label";
			current_resolution = 1;
			current_full_sizeX = sizeX;
			current_full_sizeX = sizeZ;
			current_channels = "NA";
		} else if (matches(name,"overview")) {
			current_fov = "overview";
			current_type = "overview";
			current_resolution = 1;
			current_full_sizeX = sizeX;
			current_full_sizeZ = sizeZ;
			current_channels = "NA";
		} else if (matches(name,"macro image")) {
			current_fov = "macro image";
			current_type = "macro image";
			current_resolution = 1;
			current_full_sizeX = sizeX;
			current_full_sizeZ = sizeZ;
			current_channels = "NA";
		} else if (matches(name, ".*vsi #.*")) {
			current_resolution++;
		} else if (matches(name,".*x_.*")) {
			parts = split(name,"_");
			parts2 = split(parts[2],"#");
			current_mag = parts[0];
			current_fov = parseInt(parts2[0]);
			current_channels = parts[1];
			current_resolution = 1;
			current_type = name;
			current_index++;
			current_full_sizeX = sizeX;
			current_full_sizeZ = sizeZ;			
		}
		info[m*i+0] = name;
		info[m*i+1] = current_type;
		info[m*i+2] = current_resolution;
		info[m*i+3] = current_pixel_size_x * current_full_sizeX / sizeX;
		info[m*i+4] = current_pixel_size_z * current_full_sizeZ / sizeZ;
		info[m*i+5] = current_mag;
		info[m*i+6] = current_fov;
		info[m*i+7] = "\""+current_channels+"\"";
		info[m*i+8] = sizeX;
		info[m*i+9] = sizeY;
		info[m*i+10] = sizeZ;
		info[m*i+11] = sizeC;
		if (matches(current_type,"label") || matches(current_type,"overview") || matches(current_type,"macro image"))  {
			info[m*i+12] = "0";
		} else {
			info[m*i+12] = current_index;
		}
	}
	return info;
}

function showCellSenseInfo(info) {
	/*
	 * Display information for each serie (resolution level) in the vsi file
	 * info: array as 10 fields x Nseries
	 */
	Table.create("Cell Sense");
	m = 13;
	n = info.length / m;
	for (i = 0; i < n; i++) {
		Table.set("Serie",i,          i);
		Table.set("Name",i,          info[m*i+0]);
		Table.set("Type",i,          info[m*i+1]);
		Table.set("Resolution",i,    info[m*i+2]);
		Table.set("Pixel size XY",i, info[m*i+3]);
		Table.set("Pixel size Z",i,  info[m*i+4]);
		Table.set("Magnification",i, info[m*i+5]);
		Table.set("FOV",i,           info[m*i+6]);
		Table.set("Channels",i,      info[m*i+7]);
		Table.set("Width",i,         info[m*i+8]);
		Table.set("Heigth",i,        info[m*i+9]);
		Table.set("Depth",i,         info[m*i+10]);
		Table.set("Channels",i,      info[m*i+11]);
		Table.set("Index",i,         info[m*i+12]);
	}
	Table.update;
}

function listCellSensImages(info) {
	/* Count the number of images (besides label, overview of macro image) */
	n = 0;
	m = 13;
	indices = newArray(info.length / m);
	J = 0;
	for (i = 0; i <  info.length / m; i++) {
		idx = info[m*i+12];
		if (!matches(idx, "0")) {
			newidx = true;
			for (j = 0; j < J; j++) {
				if (idx == indices[j]) {
					newidx = false;
					break;
				}
			}
			if (newidx) {
				indices[J] = idx;
				J++;				
			}			
		}
	}
	return Array.trim(indices, J);
}

function getCellSenseSerie(info, image, resolution) {
	/* Get the serie index corresponding to the given image and resolution */
	m = 13;
	for (i = 0; i < info.length / m; i++){
		if (info[m*i+12] == image && info[m*i+2] == resolution) {
			return i;
		}
	}
	return -1;
}

function getCellSenseValue(info, image, resolution, field) {
	m = 13;
	i = getCellSenseSerie(info, image, resolution);
	keys = newArray("Name","Type","Resolution","Pixel size XY","Pixel size Z","Magnification","FOV","Channels","Width","Height","Depth","Channels","Index");
	values = Array.getSequence(m);
	List.fromArrays(keys, values);
	j = List.getValue(field);
	List.clear();
	return info[m*i+j];
}

function loadCellSensVSIImage(path, info, image, resolution) {
	/* load a cell sense image */	
	serie = getCellSenseSerie(info, image, resolution);	
	width = getCellSenseValue(info, image, resolution, "Width");
	height = getCellSenseValue(info, image, resolution, "Height");
	depth = getCellSenseValue(info, image, resolution, "Depth");
	channels = getCellSenseValue(info, image, resolution, "Channels");
	title = getCellSenseValue(info, image, resolution, "Type");
	pixel_size_xy = getCellSenseValue(info, image, resolution, "Pixel size XY");
	pixel_size_z = getCellSenseValue(info, image, resolution, "Pixel size Z");
	print("Loading image serie  " + serie + " : " + title + ":  (" + width +","+height+","+channels+","+depth+")");
	newImage(title, "16-bit composite-mode", width, height, channels, depth, 1);
	id = getImageID();
	
	Ext.setSeries(serie);	
	for (i = 0; i < depth * channels; i++) {		
		Ext.getZCTCoords(i, z, c, t);				
		Ext.openImage("plane", i);
		run("Copy");
		close();
		selectImage(id);
		Stack.setPosition(c+1, z+1, t+1);
		run("Paste");		
	}
	run("Select None");

	for (channel = 1; channel <= channels; channel++) {
		Stack.setChannel(channel);
		resetMinAndMax();
	}
	setVoxelSize(pixel_size_xy, pixel_size_xy, pixel_size_z, "micron");
	return id;
}

function loadCellSensVSITile(path, info, image, resolution, x, y, tile_width, tile_height) {
	/* Open image tiles. */
	serie = getCellSenseSerie(info, image, resolution);		
	depth = getCellSenseValue(info, image, resolution, "Depth");
	channels = getCellSenseValue(info, image, resolution, "Channels");
	title = getCellSenseValue(info, image, resolution, "Type");
	pixel_size_xy = getCellSenseValue(info, image, resolution, "Pixel size XY");
	pixel_size_z = getCellSenseValue(info, image, resolution, "Pixel size Z");		
	print("Loading image serie  " + serie + " : " + title + ":  (" + tile_width +","+tile_height+","+channels+","+depth+")");
	newImage(title+"_x"+x+"_y"+y, "16-bit composite-mode", tile_width, tile_height, channels, depth, 1);
	id = getImageID();
	Ext.setSeries(serie);	
	for (i = 0; i < depth * channels; i++) {		
		Ext.getZCTCoords(i, z, c, t);	
		Ext.openSubImage("plane", i, x, y, w, h);		
		run("Copy");
		close();
		selectImage(id);
		Stack.setPosition(c+1, z+1, t+1);
		run("Paste");		
	}
	run("Select None");

	for (channel = 1; channel <= channels; channel++) {
		Stack.setChannel(channel);
		resetMinAndMax();
	}
	setVoxelSize(pixel_size_xy, pixel_size_xy, pixel_size_z, "micron");	
	return id;
}


function saveImage(opath, ipath, name, tag, format) {	
	Stack.getDimensions(width, height, channels, slices, frames);
	ofilename =  opath + File.separator + File.getNameWithoutExtension(ipath) + "_" + replace(name, ",", "_") + tag;		
	if (format=="PNG") {
		if (channels > 1) {
			run("RGB Color");
		} else {
			run("8-bit");
		}		
		saveAs("PNG", ofilename + ".png");
		if (channels > 1) {
			close();
		}
	} else if (format=="JPEG") {
		if (channels > 1) {
			run("RGB Color");
		} else {
			run("8-bit");
		}
		saveAs("JPEG", ofilename + ".jpg");
		if (channels > 1) {
			close();
		}
	} else if (format=="TIFF") {
		saveAs("TIFF", ofilename + ".tif");
	}
}

function collectGarbage() {
	for (i = 0; i < 5; i++) {
		run("Collect Garbage");
	}
}

function setLutByWavelength(w)  {
	run("Set Color By Wavelength", "wavelength="+w);
	v = getValue("rgb.foreground");
	p = newArray((v>>16)&0xff, (v>>8)&0xff, v&0xff);	
	r = newArray(256);
	g = newArray(256);
	b = newArray(256);
	for (i= 0 ; i < 256; i++) {
		r[i] = p[0]*i / 255;
		g[i] = p[1]*i / 255;
		b[i] = p[2]*i / 255;
	}
	setLut(r,g,b);
}

function applyLuts(lutstr) {
	Stack.getDimensions(width, height, channels, slices, frames);
	lut_list = split(lutstr,",");
	for (i = 0; i < lut_list.length; i++) {
		lut_list[i] = String.trim(lut_list[i]);
	}
	for (channel = 1; channel <= channels; channel++) {
		Stack.setChannel(channel);
		if (channel - 1 < lut_list.length) {			
			if (matches(lut_list[channel-1], "DAPI")) {
				setLutByWavelength(457);
			} else if (matches(lut_list[channel-1], "FITC")) {
				setLutByWavelength(512);
			} else if (matches(lut_list[channel-1], "Cy3")) {
				setLutByWavelength(569);
			} else if (matches(lut_list[channel-1], "Cy5")) {
				setLutByWavelength(670);
			} else if (matches(lut_list[channel-1], "Cy7")) {
				setLutByWavelength(756);
			} else {
				run(lut_list[channel-1]);
			}
		}
	}	
}

function exportCellSensVSIImageList(opath, ipath, info, resolution, luts, subset, projection, tag, format) {
	images = listCellSensImages(info);
	setBatchMode("hide");
	print(luts);
	for (i = 0; i < images.length; i++) {		
		id = loadCellSensVSIImage(ipath, info, images[i], resolution);
		if (luts != "None") {
			applyLuts(luts);
		}
		if (subset != "all") {
			old = getImageID();
			run("Make Subset...", subset);
			selectImage(old); close();
		}
		if (projection != "None") {		
			old = getImageID();
			run("Z Project...", "projection=["+projection+"]");
			selectImage(old); close();
		} 					
		name = getCellSenseValue(info, images[i], resolution, "Type");	
		saveImage(opath, ipath, name, tag, format);	
		close();
		collectGarbage();
	}
	setBatchMode("exit and display");
}


function processellSensVSImageByTile(path, info, image, resolution, tile_width, tile_height, luts, subset, projection, macro_name) {
	/* Process the image in tiles */	
	width = getCellSenseValue(info, image, resolution, "Width");
	height = getCellSenseValue(info, image, resolution, "Height");	
	serie = getCellSenseSerie(info, image, resolution);
	pixel_size_xy = getCellSenseValue(info, image, resolution, "Pixel size XY");
	pixel_size_z = getCellSenseValue(info, image, resolution, "Pixel size Z");
	Ext.setId(path);
	Ext.setSeries(serie);
	coords = newArray(0);
	for (y = 0; y < height; y += tile_height) {
		for (x = 0; x < width; x += tile_width) {
			w = minOf(x + tile_width, width) - x;
			h = minOf(y + tile_height, height) - y;
			id  = loadCellSensVSITile(ipath, info, image, resolution, x, y, w, h);			
			if (luts != "None") {
				applyLuts(luts);
			}
			if (subset != "all") {
				old = getImageID();
				run("Make Subset...", subset);
				selectImage(old); close();
			}
			if (projection != "None") {		
				old = getImageID();
				run("Z Project...", "projection=["+projection+"]");
				selectImage(old); close();
			} 			
			runMacro(macro_name, "x="+x+" y="+y);
		}
	}	
}

function processCellSensVSIImageList(path, info, resolution_level, tile_width, tile_height, luts, subset, projection, macroname) {
	/* Process each image in the cellsens VSI file */
	
	images = listCellSensImages(info);
	setBatchMode("hide");
	print(luts);
	for (i = 0; i < images.length; i++) {							
		processellSensVSImageByTile(path, info, images[i], resolution_level, tile_width, tile_height, luts, subset, projection, macroname);		
		collectGarbage();
	}
	setBatchMode("exit and display");	
}

start_time = getTime();
run("Bio-Formats Macro Extensions");

info = parseCellSenseSeries(ipath);

if (action == "Show information") {		
	showCellSenseInfo(info);	
} else if (action == "Export Images") {
	exportCellSensVSIImageList(opath, ipath, info, resolution_level, luts, subset, projection, tag, format);
} else if (action == "Process Images") {
	processCellSensVSIImageList(ipath, info, resolution_level, tile_size, tile_size, luts, subset, projection, macroname);
}

Ext.close();
end_time = getTime();
print("Elapsed time " + (end_time - start_time) / 1000 + " seconds.");
