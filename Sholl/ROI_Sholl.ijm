//@File (label="Input folder", value="/media/cephfs/kneish/Ctl Straitum 2") input_folder
//@File (label="Output folder", value="/media/cephfs/kneish/Ctl Straitum 2") output_folder
//@String (label="image file extension", value=".ims") imgext
//@String (label="roi file extension", value=".zip") roiext
//@Integer (label="Channel", value=4) channel
//@Float (label="Minimum distance [micron]",value=1) min_radius 
//@Float (label="Maximum distance [micron]",value=25) max_radius
//@Float (label="Step distance [micron]",value=0.5) step_radius
//@Boolean (label="Use SNT",value=true) use_snt

/*
 * Automate sholl analysis 
 * 
 * Image and ROI are saved in individual folder
 * sq
 * The macro will
 * - open the image and ROI in a folder
 * - compute a maximum intensity projection
 * - fore each ROI, segment neurites and perform a sholl analysis
 * 
 * https://imagej.net/plugins/snt/sholl#batch-processing
 * 
 * Jerome Boulanger for Karen Neish 2025
 */
 

// Find the image file

function openImageAndROIFile(folder, imgext, roiext) {
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
	// list files in folder 
	filelist = getFileList(folder);
	imgpath = 0;
	roipath = 0;
	// find image and roi files
	for (i = 0; i < lengthOf(filelist); i++) {
	    if (endsWith(filelist[i], imgext)) {
	    	imgpath = folder + File.separator + filelist[i];
	    }
	    if (endsWith(filelist[i], roiext)) {
	    	roipath = folder + File.separator + filelist[i];
	    }
	    // stop looking if the both are found
	    if (imgpath != 0 && roipath !=0) {
	    	break;
	    }
	}
	// load image and roi files
	print("Opening image file:\n" + imgpath);
	run("Bio-Formats Importer", "open=["+imgpath+"] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_1");
	id = getImageID();
	print("Opening ROI file:\n" + roipath);
	open(roipath);
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
	for (i = 0; i < x.length; i++) {
		x[i] *= factor;
		y[i] *= factor;
	}
	Roi.setPolygonSplineAnchors(x, y);
	run("Select None");
	run("Restore Selection");
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
	getDimensions(width0, height, channels, slices, frames);
	for (iter = 0; iter < 6; iter++) {
		selectImage(id);
		zoom(1.25);
		run("FeatureJ Hessian", "smallest smoothing=0.25");
		run("Macro...", "code=v=-v*(v<0)");
		nid = getImageID();
		imageCalculator("Add", id, nid);
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
}

function segmentNeurite(id) {
	/*
	 * Segment neurites in the image
	 * 
	 * @param (id): image id 
	 */
	run("FeatureJ Hessian", "smallest smoothing=0.25");
	feat = getImageID();
	setLaplacianThreshold(feat, 0.01);
	run("Skeletonize");
	return getImageID();
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
	N = Math.ceil(2 * PI * radius);
	array = newArray(3 * N);
	for (i = 0 ; i < N; i++) {
		array[i+n] = (center[0] + radius * cos(i * 2 * PI / N)) / dx;
		array[i+2*n] =  (center[1] + radius * sin(i * 2 * PI / N)) / dy;
		array[i] = getPixel(array[i+n], array[i+2*n]);
	}
	return array;
}

function linspace(min_radius, max_radius, step_radius) {
	n = (max_radius - min_radius) / step_radius + 1;
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

function getSomaCenter(id) {
	/*
	 * compute the center in micron
	 */
	selectImage(id);
	run("Duplicate...","title=tmp");
	getStatistics(area, mean, min, max, std, histogram);
	setThreshold(mean+std,max);
	run("Convert to Mask");
	run("Minimum...", "radius=15");
	run("Maximum...", "radius=15");
	run("Create Selection");
	center_x = getValue("X");
	center_y = getValue("Y");
	close();
	return newArray(center_x, center_y);
}

function processROI(id, i, min_radius, max_radius, step_radius, tbl, output_folder) {
	/*
	 * Process a single ROI
	 * 
	 * Store results in a table and export crop to a folder
	 * 
	 * @param (int): image id
	 * @param (int): roi index
	 */
	setBatchMode("hide");
	// select the original image
	selectImage(id);
	getPixelSize(unit, dx, dy);
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
	center = getSomaCenter(id);
	print("   x:"+center[0]+"um, y:"+center[1]+"um");
	//setBatchMode("exit and display");
	distances = linspace(min_radius, max_radius, step_radius);
	if (!isOpen(tbl)) {
		Table.create(tbl);
		Table.setColumn("distance", distances);
	}
	// sholl analysis
	if (use_snt) {
		selectImage(skel);
		makePoint(center[0]/dx, center[1]/dy);
		run("Legacy: Sholl Analysis (From Image)...", "starting="+min_radius+" ending="+max_radius+" radius_step="+step_radius+" #_samples=1 integration=Mean enclosing=1 #_primary=[] infer fit linear polynomial=[Best fitting degree] linear normalizer=Area directory=["+output_folder+"]");
		selectWindow("skeleton_Sholl-Profiles");
		Table.save(output_folder + File.separator + name + "-sholl-profiles.csv");
		number_intersection = Table.getColumn("Inters.");
		selectWindow("Sholl profile (Linear) for skeleton");
		saveAs("JPG", output_folder + File.separator + name + "-linear.jpg");
		close();
	} else {
		print("- Sholl analysis");
		selectImage(skel);
		distances = linspace(min_radius, max_radius, step_radius);
		number_intersection = shollAnalysis(center, distances);
	}

	selectWindow(tbl);
	Table.setColumn(name, number_intersection);
	
	// Save results
	print("- Create composite image");
	selectImage(enhanced);
	run("8-bit");
	getDimensions(width, height, channels, slices, frames);
	selectImage(crop);
	run("Size...", "width="+width+" height="+height+" depth=1 average interpolation=None");
	run("Merge Channels...", "c1=crop c2=enhanced c3=skeleton create");
	merge = getImageID();
	run("Restore Selection");
	Overlay.addSelection
	
	makePoint(center[0]/dx, center[1]/dy);
	Overlay.addSelection
	ofile = output_folder + File.separator + name + ".tif";
	print("- Saving cropped image");
	print(ofile);
	saveAs("TIFF", ofile);
	selectImage(merge); close();
	run("Collect Garbage");

	setBatchMode("exit and display");
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

function main(input_folder, output_folder, imgext, roiext, channel, min_radius, max_radius, step_radius,tbl) {
	start_time = getTime();
	// close all image windows
	run("Close All");
	// close the roi manager
	if (isOpen("ROI Manager")) {
		selectWindow("ROI Manager");
		run("Close");
	}
	// closr the table
	if (isOpen("Sholl Analysis")) {
		selectWindow("Sholl Analysis");
		run("Close");
	}
	// close sholl result
	if (isOpen("Sholl Results")) {
		selectWindow("Sholl Results");
		run("Close");
	}
	// open the image
	id = openImageAndROIFile(input_folder, imgext, roiext);
	// z projection and channel selection
	mip = projectAndSelectChannel(id, channel);
	// process the ROIs (crop+preprocess+Sholl analysis+export results)
	processAllRoi(mip, min_radius, max_radius, step_radius,tbl, output_folder);
	// export the table
	selectWindow("Sholl Analysis");
	Table.save(output_folder + File.separator + "sholl-analysis.csv");
	if (use_snt) {
		selectWindow("Sholl Results");
		Table.save(output_folder + File.separator + "sholl-results.csv");
	}
	end_time = getTime();
	print("Finished in " + (end_time - start_time)/1000 + " seconds.");

}

main(input_folder, output_folder, imgext, roiext, channel, min_radius, max_radius, step_radius,"Sholl Analysis");

