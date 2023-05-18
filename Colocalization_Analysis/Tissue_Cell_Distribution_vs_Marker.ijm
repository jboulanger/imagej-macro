#@File   (label="Input file (VSI)", description="Select a cellsens (VSI file)") path
#@Integer(label="Image index",                     value=1,    description="Image index") image
#@Integer(label="Cell channel",                    value=2,    description="Channel with cell marker") channel1
#@Integer(label="Resolution for cell detection",   value=3,    description="Resolution level in the VSI file") resolution1
#@Float  (label="Threshold for cell detection",    value=5000, description="Intensity threshold") threshold1
#@Integer(label="Tile size for cell detection",    value=2000, description="Tile size in pixels") tile_size
#@Integer(label="Reference marker channel",        value=4,    description="Channel for the marker") channel2
#@Integer(label="Resolution for reference marker", value=6,    description="Resolution level in the VSI file.") resolution2
#@String(label="Action", choices={"Process","Show Info","Display Cell Channel","Display Marker Channel"}) action

/*
 * Analysis of spatial organization of cell with respect to a reference marker
 *
 * Cells and marker are in different channels and segmented in different 
 * resolution from the VSI image.
 *  
 *  Open the macro in the script editor and press "run"
 *
 *  Require the IJPB-plugins
 * 
 * Jerome for Jane 2023
 */

run("Close All");
run("Bio-Formats Macro Extensions");

// convert the channel indices to bioformat extension convention 
channel1--;
channel2--;

info = parseCellSenseSeries(path);
if (matches(action, "Show Info")) {
	showCellSenseInfo(info);
	exit();
} else if (matches(action, "Display Cell Channel")) {
	showCenterCrop(path, info, image, resolution1, channel1, tile_size, tile_size);
	exit();
} else if (matches(action, "Display Marker Channel")) {
	showCenterCrop(path, info, image, resolution1, channel2, tile_size, tile_size);
	exit();
}
coords = detectBrightCells(path, info, image, resolution1, channel1, threshold1, tile_size, tile_size);
mask = segmentRegion(path, info, image, resolution2, channel2);
edm = distanceToRegion(info, image, resolution2, mask, coords);
selectImage(mask); close();
selectImage(edm);
statisticalRandomAnalysis(edm);
selectImage(edm); close();
createVisualization(path, info, image, resolution2);


function createVisualization(path, info, image, resolution) {
	/*
	 * Load the image at low res and add the ROIs as overlays
	 */
	Ext.setId(path);
	serie = getCellSenseSerie(info, image, resolution);
	Ext.setSeries(serie);
	run("Bio-Formats Importer", "open=["+path+"] autoscale color_mode=Composite rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_"+(serie+1));
	n = roiManager("count");
	for (i = 0; i <	n; i++) {
		roiManager("select", i);
		run("Add Selection...");
	}
}

function getArrayHistogram(a,nbins,min,max) {
	/*
	 * Histogram of the values of the array a
	 * 
	 * a     : array
	 * nbins : number of bins
	 * min   : minimum range value
	 * max   : maximum range value
	 */
	h = newArray(nbins);
	for (i = 0; i < a.length; i++) {
		v = a[i];
		if (v >= min && v < max) {
			k = round( (v - min) / (max - min) *  (nbins - 1));
			h[k] = h[k] + 1.0 / a.length;
		}
	}
	return h;
}

function statisticalRandomAnalysis(edm) {
	/*
	 * Statistical analysis of distance distribution
	 * 
	 * edm: ID of the eucidean distance map
	 * 
	 */
	selectImage(edm);
	roiManager("select", 0);
	x0 = getValue("BX");
	y0 = getValue("BY");
	width = getValue("Width");
	height = getValue("Height");
	max = getValue("Max");
	nbins = 30;
	getHistogram(radius, counts, nbins, 0, max);
	s = 0;
	for (i = 0; i < counts.length; i++) {
		s += counts[i];
	}
	reference = counts;
	for (i = 0; i < reference.length; i++) {
		reference[i] = counts[i] / s;
	}
	selectWindow("Positions");
	dist = Table.getColumn("Distance [um]");
	measure = getArrayHistogram(dist, nbins, 0, max);
	Array.show("Histograms", radius, reference, measure);

	//	randomization test
	N = dist.length;// number of points
	K = 100; // number of trial
	sim = newArray(N*K);
	for (k = 0; k < K; k++) { // number of trial
		selectImage(edm);
		roiManager("select", 0);		
		for (i = 0; i < N; i++) {
			// generate a random point in the ROI 0
			for (j = 0; j < 1000; j++) {
				x = x0 + width * random;
				y = y0 + height * random;
				if (Roi.contains(x,y)) {
					sim[i+N*k] = getPixel(x, y);
					break;
				}
			}
		}
		
	}
	print("Creating figure");
	Plot.create("Distance analysis", "Distance [um]", "Frequency");
	
	Plot.setColor("#AAAAAA");
	for (k = 0; k < K; k++) {
		simk = Array.slice(sim,k*N,(k+1)*N-1);
		hk = getArrayHistogram(simk, nbins, 0, max);
		Plot.add("line", radius, hk);
	}
	
	Plot.setColor("blue");
	Plot.add("line", radius, reference);
	Plot.setColor("red");
	Plot.add("line", radius, measure);
	Plot.setLimitsToFit();
	Plot.setLegend("Reference\nMeasure", "top-right");
	Plot.update();
}

function distanceToRegion(info, image, resolution, maskid, coords) {
	/*
	 * Record the distance from the mask of each coordinates in coords.
	 * Add Rois
	 */
	dx = getCellSenseValue(info, image, resolution, "Pixel size");
	selectImage(maskid);
	run("Select None");
	run("Invert");
	run("Chamfer Distance Map", "distances=[Quasi-Euclidean (1,1.41)] output=[32 bits] normalize");
	run("Multiply...","value="+dx);
	distance_id = getImageID();
	d = 2; // size of the marker
	k = 0;
	Table.create("Positions");
	for (i = 0; i < coords.length/2; i++) {
		roiManager("select", 0);
		x = coords[2*i] / dx;
		y = coords[2*i+1] / dx;
		if (Roi.contains(x, y)) {
			makeOval(x-d/2, y-d/2, d, d);
			Roi.setStrokeColor("blue");
			roiManager("add");
			Table.set("X [um]", k, x);
			Table.set("Y [um]", k, y);
			Table.set("Distance [um]", k, getValue("Mean"));
			k++;
		}
	}
	//roiManager("Show All without labels");
	roiManager("Show None");
	return distance_id;
}

function segmentRegion(path, info, image, resolution, channel) {
	/*
	 * 	 Segment a reference large region
	 * 	 
	 * 	 path: path to the VSI file
	 * 	 info: 
	 * 	 return a mask and populate the roi manager
	 * 	 
	 */
	if (isOpen("ROI Manager")) {selectWindow("ROI Manager"); run("Close"); }
	serie = getCellSenseSerie(info, image, resolution);
	dx = getCellSenseValue(info, image, resolution, "Pixel size");
	Ext.setId(path);
	Ext.setSeries(serie);
	Ext.openImage("low res", channel);
	rename("mask");
	id = getImageID();

	// create a mask
	run("Duplicate...","duplicate");
	id0 = getImageID();
	run("Select None");
	run("Square Root");
	for (iter = 0; iter < 5; iter++) {
		run("Median...", "radius=5");
		run("Minimum...", "radius=5");
		run("Maximum...", "radius=5");
	}
	setAutoThreshold("Otsu dark");
	//setThreshold(getValue("Mean"), getValue("Max")+1);
	run("Convert to Mask");
	run("Create Selection");
	roiManager("add");
	selectImage(id0);close();

	// segment bright areas and compute distance map
	selectImage(id);
	run("Select None");
	run("Subtract Background...", "rolling=10");
	roiManager("select", 0);
	run("Clear Outside");
	run("Select None");
	setAutoThreshold("MaxEntropy dark");
	//setThreshold(threshold, 65635);
	run("Convert to Mask");
	run("Analyze Particles...", "size=0-Infinity add");
	n = roiManager("count");
	for (i = 1; i < n; i++) {
		roiManager("select", i);
		Roi.setStrokeColor("white");
		roiManager("update");
	}

	return getImageID();
}

function showCenterCrop(path, info, image, resolution, channel, tile_width, tile_height) {
	width = getCellSenseValue(info, image, resolution, "X");
	height = getCellSenseValue(info, image, resolution, "Y");
	print("Image size " + width + " x " + height);
	serie = getCellSenseSerie(info, image, resolution);
	dx = getCellSenseValue(info, image, resolution, "Pixel size");
	Ext.setId(path);
	Ext.setSeries(serie);
	x = maxOf(0, round(width/2 - tile_width / 2));
	y = maxOf(0, round(height/2 - tile_height / 2));
	w = minOf(x + tile_width, width) - x;
	h = minOf(y + tile_height, height) - y;
	Ext.openSubImage("test", channel, x, y, w, h);
}

function detectBrightCells(path, info, image, resolution, channel, threshold, tile_width, tile_height) {
	/*
	 * Detect bright cells in given file and return their coordinates as an array
	 */
	setBatchMode(true);
	width = getCellSenseValue(info, image, resolution, "X");
	height = getCellSenseValue(info, image, resolution, "Y");
	print("Image size " + width + " x " + height);
	serie = getCellSenseSerie(info, image, resolution);
	dx = getCellSenseValue(info, image, resolution, "Pixel size");
	Ext.setId(path);
	Ext.setSeries(serie);
	coords = newArray(0);
	for (y = 0; y < height; y += tile_height) {
		for (x = 0; x < width; x += tile_width) {
			w = minOf(x + tile_width, width) - x;
			h = minOf(y + tile_height, height) - y;
			print("Tile : " + w + "x" + h);
			ncoords = processTile(x,y,w,h,threshold, dx);
			coords = Array.concat(coords, ncoords);
		}
	}
	setBatchMode(false);
	return coords;
}


function processTile(x,y,w,h, threshold, dx) {
	/*
	 * Open and detect bright cells in tiles.
	 */
	 
	// open the tile
	Ext.openSubImage("test", channel, x, y, w, h);
	id = getImageID();
	// process the image
	//run("Subtract Background...", "rolling=100");
	setThreshold(threshold, getValue("Max") + 1);
	if (isOpen("ROI Manager")) { selectWindow("ROI Manager"); run("Close"); }
	run("Select None");
	run("Analyze Particles...", "size=0-Infinity pixel clear add");
	// Store the coordinates of the ROI in an array
	coords = getROIsCenterCoordinates(x, y, dx);
	// close the tile
	selectImage(id); close();
	run("Collect Garbage");

	return coords;
}

function getROIsCenterCoordinates(x0,y0,dx) {
	/*
	 * Store the coordinated of the centers in a array
	 */
	n = roiManager("count");
	coords = newArray(2*n);
	getPixelSize(unit, pixelWidth, pixelHeight);
	for (i = 0; i < n; i++) {
		roiManager("select", i);
		coords[2*i] = x0 + dx * getValue("X") / pixelWidth;
		coords[2*i+1] = y0 + dx * getValue("Y") / pixelHeight;
	}
	return coords;
}

function parseCellSenseSeries(path) {
	/*
	 * Parse file and return an info array with
	 * 0: name,
	 * 1: type label/overview/image/macroimage
	 * 2: resolution level (1-6)
	 * 3: pixel size
	 * 3: objective mag
	 * 4: field of view
	 * 5:
	 */
	m = 12; // number of fields
	Ext.setId(path);
	Ext.getSeriesCount(n);
	info = newArray(m*n);
	current_type = "";
	current_resolution = 1;
	current_fov = "";
	current_mag = "";
	current_index = 0;
	current_full_sizeX = 1;
	for (i = 0; i < n; i++) {
		Ext.setSeries(i);
		Ext.getSeriesName(name);
		Ext.getSizeX(sizeX);
		Ext.getSizeY(sizeY);
		Ext.getSizeZ(sizeZ);
		Ext.getSizeC(sizeC);
		Ext.getSizeT(sizeT);
		Ext.getPixelsPhysicalSizeX(dx);
		current_pixel_size = dx;
		if (matches(name,"label")) {
			current_fov = "label";
			current_type = "label";
			current_resolution = 1;
			current_full_sizeX = sizeX;
			current_channels = "NA";
		} else if (matches(name,"overview")) {
			current_fov = "overview";
			current_type = "overview";
			current_resolution = 1;
			current_full_sizeX = sizeX;
			current_channels = "NA";
		} else if (matches(name,"macro image")) {
			current_fov = "macro image";
			current_type = "macro image";
			current_resolution = 1;
			current_full_sizeX = sizeX;
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
		}
		info[m*i] = name;
		info[m*i+1] = current_type;
		info[m*i+2] = current_resolution;
		info[m*i+3] = current_pixel_size * current_full_sizeX / sizeX;
		info[m*i+4] = current_mag;
		info[m*i+5] = current_fov;
		info[m*i+6] = "\""+current_channels+"\"";
		info[m*i+7] = sizeX;
		info[m*i+8] = sizeY;
		info[m*i+9] = sizeZ;
		info[m*i+10] = sizeC;
		if (matches(current_type,"label") || matches(current_type,"overview") || matches(current_type,"macro image"))  {
			info[m*i+11] = "0";
		} else {
			info[m*i+11] = current_index;
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
	m = 12;
	n = info.length / m;
	for (i = 0; i < n; i++) {
		Table.set("Serie",i,          i);
		Table.set("Name",i,          info[m*i]);
		Table.set("Type",i,          info[m*i+1]);
		Table.set("Resolution",i,    info[m*i+2]);
		Table.set("Pixel size",i,    info[m*i+3]);
		Table.set("Magnification",i, info[m*i+4]);
		Table.set("FOV",i,           info[m*i+5]);
		Table.set("Channels",i,      info[m*i+6]);
		Table.set("X",i,             info[m*i+7]);
		Table.set("Y",i,             info[m*i+8]);
		Table.set("Z",i,             info[m*i+9]);
		Table.set("C",i,             info[m*i+10]);
		Table.set("Index",i,         info[m*i+11]);
	}
	Table.update;
}

function numberOfCellSenseImages(info) {
	/*
	 * Count the number of images (besides label, overview of macro image)
	 */
	n = 0;
	m = 12;
	for (i = 0; i <  info.length / m; i++) {
		if (!matches(info[m*i+11], "0")) {
			n++;
		}
	}
	return n;
}

function getCellSenseSerie(info, image, resolution) {
	/*
	 * Get the serie index corresponding to the given image and resolution
	 */
	m = 12;
	for (i = 0; i < info.length / m; i++){
		if (info[m*i+11] == image && info[m*i+2] == resolution) {
			return i;
		}
	}
	return -1;
}

function getCellSenseValue(info, image, resolution, field) {
	m = 12;
	i = getCellSenseSerie(info, image, resolution);
	keys = newArray("Name","Type","Resolution","Pixel size","Magnification","FOV","Channels","X","Y","Z","C","Index");
	values = Array.getSequence(m);
	List.fromArrays(keys, values);
	j = List.getValue(field);
	List.clear();
	return info[m*i+j];
}
