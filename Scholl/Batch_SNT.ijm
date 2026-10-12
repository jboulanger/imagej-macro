//@File (label="Input mage") input_path
//@File (label="Output folder", style="directory") output_folder
//@String (label="roi file extension", value=".zip") roiext
//@Integer (label="Channel", value=4) channel
//@Float (label="Minimum distance [micron]",value=1) min_radius 
//@Float (label="Maximum distance [micron]",value=25) max_radius
//@Float (label="Step distance [micron]",value=0.5) step_radius

/*
 * Automate sholl analysis for multiple ROI stored in a ZIP file
 * 
 * Image and ROI are saved in individual folder
 * Results are saved in an subfolder.
 * 
 * The macro will
 * - open the image and corresponding ROI file 
 * - compute a maximum intensity projection and extract the channel with neurites
 * - for each ROI, segment neurites and perform a sholl analysis
 * 
 * Jerome Boulanger for Karen Neish 2025, Nicole Morissey 2026
 */
 

// Find the image file

function openImageAndROIFile(image_path, roiext) {
	/*
	 * Inspect the folder to find an image file with imgext and roi file with roiext
	 * 
	 * Open the first serie of the image file using bioformat
	 * 
	 * 	@param folder: path to the folder
	 * 	@param imgext: extension of the image file
	 * 	@param roiext: extension of the roi file
	 * 	@return id of the image window
	 * 	
	 */
	 
	// Open the image
	print("Opening image file:\n" + image_path);
	run("Bio-Formats Importer", "open=["+image_path+"] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_1");
	id = getImageID();
	
	// Locates and open the corresponding ROI
	folder = File.getDirectory(image_path);
	roi_path = folder + File.separator + File.getNameWithoutExtension(image_path) + roiext;
	
	print("Opening ROI file:\n" + roi_path);
	open(roi_path);
	return id;
}

function projectAndSelectChannel(id, channel) {
	/* Maximum intensity projection and channel selection
	 *  
	 * @param id: window id
	 * @param channel: index of the channel
	 * @return : id of the window with mip+channel selection
	 */
	setBatchMode(true);
	selectImage(id);
	print("- Maximum intensity projection");
	run("Z Project...", "projection=[Max Intensity]");
	mipid = getImageID();
	print("- Extracting channel "+ channel);
	run("Duplicate...", "duplicate channels="+channel);
	chid = getImageID();
	selectImage(mipid); close();
	selectImage(id); close();
	selectImage(chid);
	setBatchMode(false);
	return chid;
}

function image2Array() {
	/*
	 * Load pixel values in an array
	 * 
	 * @return (array): pixels intensities
	 */
	getDimensions(width, height, channels, slices, frames);
	a = newArray(width * height);
	for (y = 0; y < height; y++) {
		for (x = 0;x < width; x++) {
			a[x+width*y] = getPixel(x, y);
		}
	}
	return a;
}

function subtractConstantBackground() {
	/*
	 * Subtract a constant background
	 * 
	 * The background is estimated as the 10th percentile.
	 */
	run("Select None");
	pixels = image2Array();
	Array.sort(pixels);
	bg = pixels[round(pixels.length / 10)];
	run("Subtract...", "value="+bg);
	run("Restore Selection");
}

function deleteOffSection() {
	/*
	 * Erase the content of the image outside the selection
	 */
	run("Make Inverse");
	run("Cut");
	run("Make Inverse");
}

function zoom(factor) {
	/*
	 * Scale the image by a factor in place
	 * 
	 * @param (float): zoom factor
	 */
	getDimensions(width, height, channels, slices, frames);
	run("Size...", "width="+factor*width+" height="+factor*height+" depth=1 constrain average interpolation=Bilinear");
}

function zoomROI(factor) {
	/*
	 * Zoom a ROI 
	 */
	Roi.getSplineAnchors(x, y);
	x = Array.resample(x, 6 * x.length);
	y = Array.resample(y, 6 * y.length);
	
	for (i = 0; i < x.length; i++) {
		x[i] *= factor;
		y[i] *= factor;
	}
	Roi.setPolygonSplineAnchors(x, y);
	run("Select None");
	run("Restore Selection");
}

function stepSharp(id) {
	selectImage(id);
	id = getImageID();
	run("Duplicate...", " ");
	u = getImageID();
	for (i=0;i<4;i++) {
		selectImage(u);
		run("Duplicate...", "title=hu");
		hu = getImageID();
		run("Gaussian Blur...", "sigma=2");
		run("Macro...", "code=v=1/v");
		imageCalculator("Multiply stack", hu, id);
		run("Gaussian Blur...", "sigma=1");
		imageCalculator("Multiply stack", u, hu);
		
		selectImage(u);
		run("FeatureJ Hessian", "smallest smoothing=0.25");
		hess = getImageID();
		run("Macro...", "code=v=-0.1*v*((v-10)<0)");
		imageCalculator("Add stack", u, hess);
		selectImage(hess); close();
		selectImage(hu);close();
	}
	selectImage(u);
	
}

function enhanceNeurites(id) {
	/*
	 * Enhance neurites in place
	 * 
	 * @param id: image id
	 * 
	 */
	selectImage(id);
	run("32-bit");
	run("Select None");
	//run("Subtract Background...", "rolling=20");
	//run("Square Root");
	getDimensions(width0, height, channels, slices, frames);
	for (iter = 0; iter < 4; iter++) {
		selectImage(id);
		zoom(1.25);
		run("FeatureJ Hessian", "smallest smoothing=0.25");
		getStatistics(area, mean, min, max, std, histogram);
		run("Macro...", "code=v=-10*v*((v+"+std+")<0)");
		//run("Macro...", "code=v=-5*v");
		nid = getImageID();
		imageCalculator("Add", id, nid);
		//setBatchMode("exit and display");exit();
		selectImage(nid); close();
	}
	getDimensions(width1, height, channels, slices, frames);
	run("Restore Selection");
	zoomROI(width1 / width0);
}

function setLaplacianThreshold(id, p) {
	/*
	 * Threshold the image assuming a Laplacian distribution
	 * 
	 * @param (int): image id
	 * @param (float): p value
	 * 
	 */
	values = image2Array();
	Array.getStatistics(values, min, max, mean, stdDev);
	for (i = 0; i < values.length; i++) {
		m += abs(values[i]-mean);
	}
	scale = m / values.length;
	threshold = mean - scale * log( 2 - 2 * p);
	setThreshold(min, threshold);
	run("Convert to Mask");
	run("Skeletonize");
	run("Restore Selection");
	run("Enlarge...", "enlarge=-10 pixel");
	setBackgroundColor(0, 0, 0);
	run("Clear Outside");
}

function segmentNeurite(id) {
	/*
	 * Segment neurites in the image
	 * 
	 * @param (id): image id 
	 */
	run("FeatureJ Hessian", "smallest smoothing=0.5");
	feat = getImageID();
	setLaplacianThreshold(feat, 0.01);
	run("Skeletonize");
	skel = getImageID();
	run("Connected Components Labeling", "connectivity=8 type=[16 bits]");
	lbl = getImageID();
	run("Keep Largest Label");
	arbor = getImageID();
	selectImage(skel); 
	close();
	selectImage(lbl);
	close();
	return arbor;
}

function getCircleProfile(center, radius) {
	/*
	 * Get the circle profile value,x,y
	 * 
	 * @param center (array): center of the circle [x,y]
	 * @param radius (float): radius of the circle
	 * @return array: intensity on the circle at every pixels
	 */
	getPixelSize(unit, dx, dy);
	N = Math.ceil(2 * PI * radius/dx);
	array = newArray(3 * N);
	for (i = 0 ; i < N; i++) {
		array[i+n] = (center[0] + radius * cos(i * 2 * PI / N)) / dx;
		array[i+2*n] =  (center[1] + radius * sin(i * 2 * PI / N)) / dy;
		array[i] = getPixel(array[i+n], array[i+2*n]);
	}
	return array;
}

function linspace(min_radius, max_radius, step_radius) {
	n = (max_radius - min_radius) / step_radius;
	a = newArray(n);
	for (i = 0; i < n; i++) {
		a[i] = min_radius + i * step_radius;
	}
	return a;
}

function shollAnalysis(center, distances) {
	/* 
	 *  Sholl analysis
	 *  
	 *  Count the number of intersection of the skeletonized image
	 *  with circles.
	 *  
	 *  @param center (array): center in micron
	 *  @param distance (array): distances in micron
	 *  @return (array): number of intersections
	 */
	
	n = (max_radius - min_radius) / step_radius;
	number_intersection = newArray(distances.length);
	for (i = 0; i < n; i++) {
		profile = getCircleProfile(center,distances[i]);
		N = profile.length / 3;
		peaks = Array.findMaxima(Array.slice(profile, 0, N), 0.5);
		number_intersection[i] = peaks.length;
	}
	
	return number_intersection;
}

function drawShollCircles(id, center, distances, roi, x0, y0) {
	/* 
	 *  Draw Sholl analysis cicles
	 *  
	 *  
	 *  @param center (array): center in micron
	 *  @param distance (array): distances in micron
	 *  @param roi (int) : index of the roi
	 *  @param x0 (int) : x-corner of the selection in the destination
	 *  @param y0 (int) : y-corner of the selection in the destination
	 */
	 selectImage(id);
	 getPixelSize(unit, dx, dy);
	 n = (max_radius - min_radius) / step_radius;
	 for (i = 0; i < n; i++) {
	 	d = 2*distances[i]/dx;
	 	print(distances[i], d, dx);
	 	Overlay.drawEllipse(x0+center[0]/dx-d/2, y0+center[1]/dy-d/2, d, d);
		Overlay.setPosition(0, 0, roi+1);
		Overlay.show();
	 }
}

function getSomaCenter(id) {
	/*
	 * Compute the center in micron using the brightest center
	 * 
	 * @param id:  image id of the skeleton image
	 * 
	 */
	selectImage(id);
	run("Duplicate...","title=tmp");
	run("32-bit");
	run("Gaussian Blur...", "sigma=5 scaled");
	getStatistics(area, mean, min, max, std, histogram);
	setThreshold(max, max);
	run("Create Selection");
	center_x = getValue("X");
	center_y = getValue("Y");
	getPixelSize(unit, dx, dy);
	makePoint(center_x/dx, center_y/dy);
	close();
	return newArray(center_x, center_y);
}

function getSomaCenter2(id) {
	/*
	 * Get the soma center as the point the furthest from the end points.
	 * 
	 * @param id (int)  the id of the skeleton image
	 */
	selectImage(id);
	id = getImageID();
	name = getTitle();
	getVoxelSize(dx, dy, dz, unit);
	run("Duplicate...","title=marker");
	run("32-bit");
	tmp = getImageID();
	run("Convolve...", "text1=[1 1 1 1 1\n1 1 1 1 1\n1 1 25 1 1\n1 1 1 1 1\n1 1 1 1 1\n]");
	setThreshold(7140,7140);
	setOption("BlackBackground", true);
	run("Convert to Mask");
	run("Geodesic Distance Map", "marker=[marker] mask=["+name+"] distances=[Quasi-Euclidean (1,1.41)] output=[32 bits] normalize");
	gdm = getImageID();
	getStatistics(area, mean, min, max, std, histogram);
	setVoxelSize(dx, dy, dz, unit);
	setThreshold(max-1, max);
	run("Convert to Mask");
	setThreshold(128,255);
	run("Create Selection");
	center_x = getValue("X");
	center_y = getValue("Y");
	selectImage(gdm); close();
	selectImage(tmp); close();
	selectImage(id);
	makePoint(center_x/dx, center_y/dy);
	return newArray(center_x, center_y);
}

function processROI(id, i, min_radius, max_radius, step_radius, tbl, dest) {
	/*
	 * Process a single ROI
	 * 
	 * Store results in a table and export crop to a folder
	 * 
	 * @param (int): image id
	 * @param (int): roi index
	 */
	setBatchMode(true);
	// select the original image
	selectImage(id);
	getPixelSize(unit, dx0, dy0);
	imgname = getTitle();
	// crop a the ROI
	roiManager("select", i);
	name = Roi.getName;
	print("- Cropping ROI " + name);
	run("Duplicate...","title=crop");
	crop = getImageID();
	
	// enhance
	run("Duplicate...","title=enhanced");
	print("- Enhance crop ");
	subtractConstantBackground();
	deleteOffSection();
	enhanced = getImageID();
	enhanceNeurites(enhanced);
	
	// segment
	print("- Segment and skeletonize");
	skel = segmentNeurite(enhanced);
	rename("skeleton");
	
	// get the center of the soma
	print("- Computing the center of the soma");
	center = getSomaCenter(skel);
	// center = getSomaCenter2(skel);
	print("   x:"+center[0]+"um, y:"+center[1]+"um");
	print("   x:"+center[0]/dx0+"px, y:"+center[1]/dy0+"px");
	getDimensions(width, height, channels, slices, frames);
	print(width/2,height/2);
	
	// sholl analysis
	print("- Sholl analysis");
	selectImage(skel);
	distances = linspace(min_radius, max_radius, step_radius);
	number_intersection = shollAnalysis(center, distances);
	
	if (!isOpen(tbl)) {
		Table.create(tbl);
		Table.setColumn("distance", distances);
	}
	selectWindow(tbl);
	Table.setColumn(name, number_intersection);
	
	// Create a representation
	print("- Create composite image");
	selectImage(enhanced);
	getPixelSize(unit, dx, dy);
	run("8-bit");
	run("Select None");
	run("Copy"); selectImage(dest);Stack.setFrame(i+1);Stack.setChannel(1);run("Paste"); 
	
	Roi.getBounds(x0, y0, width, height);
	
	selectImage(enhanced);
	getVoxelSize(dx, dy, dz, unit);
	getDimensions(width, height, channels, slices, frames);
	selectImage(crop);
	run("Size...", "width="+width+" height="+height+" depth=1 average interpolation=None");
	run("Select None");
	run("Copy");selectImage(dest);Stack.setFrame(i+1);Stack.setChannel(2);run("Paste");
	
	selectImage(dest);
	setVoxelSize(dx, dy, dz, unit);
	
	selectImage(skel);
	run("Select None");
	run("Copy");selectImage(dest);Stack.setFrame(i+1);Stack.setChannel(3);run("Paste");
	run("Yellow");
	run("Select None");
	drawShollCircles(dest, center, distances, i, x0, y0);
	selectImage(crop); close();
	selectImage(enhanced); close();
	run("Collect Garbage");
	setBatchMode(false);
	selectImage(id);
	return id;
}

function processAllRoi(id, min_radius, max_radius, step_radius, tbl, output_folder) {
	/*
	 * Process all ROI
	 * 
	 * @param (int): image id
	 */
	n = roiManager("count");
	print("- Processing "+n+" ROI...");
	for (i = 0; i < n; i++) {
		processROI(id, i, min_radius, max_radius, step_radius, tbl, output_folder);
	}
}

function createStack() {
	/*
	 * Create a destinatation stack with the max size of the neurones
	 * 
	 */
	nroi = roiManager("count");
	wmax = 0;
	hmax = 0;
	for (i = 0; i < nroi; i++) {
		roiManager("Select",i);
		Roi.getBounds(x, y, width, height)
		wmax = maxOf(wmax, width);
		hmax = maxOf(hmax, height);
	}
	newImage("visu", "8-bit composite-mode black", 3*wmax, 3*hmax, 3, 1, nroi);
	id = getImageID();
	return id;
}


function main(input_path, output_folder, roiext, channel, min_radius, max_radius, step_radius,tbl) {
	start_time = getTime();
	// close all image windows
	run("Close All");
	// close the roi manager
	if (isOpen("ROI Manager")) {
		selectWindow("ROI Manager");
		run("Close");
	}
	// close the table
	if (isOpen("Sholl Analysis")) {
		selectWindow("Sholl Analysis");
		run("Close");
	}
	
	// open the image
	id = openImageAndROIFile(input_path, roiext);
	
	// create empty destination stack
	dest = createStack(); 
	
	// z projection and channel selection
	mip = projectAndSelectChannel(id, channel);
	
	// process the ROIs (crop+preprocess+Sholl analysis+export results)
	processAllRoi(mip, min_radius, max_radius, step_radius, tbl, dest);
	
	// Save the destination stack
	name = File.getNameWithoutExtension(input_path);
	ofile = output_folder + File.separator + name + ".tif";
	print("- Saving cropped image");
	print(ofile);
	selectImage(dest);
	saveAs("TIFF", ofile);
	
	// Save the table
	selectWindow("Sholl Analysis");
	Table.save(output_folder + File.separator + name+"-Sholl Analysis.csv");
	
	end_time = getTime();
	print("Finished in " + (end_time - start_time)/1000 + " seconds.");
}

main(input_path, output_folder, roiext, channel, min_radius, max_radius, step_radius,"Sholl Analysis");

