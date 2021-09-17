#@File    (label="Input File",description="Either a csv file or an image file") filename
#@String  (label="Group",description="Experimental condition", value="control") condition
#@String  (label="Channels name", value="DAPI,GFP") channel_names_str
#@Integer (label="ROI channel", value=1) roi_channel
#@Boolean (label="Segment Nuclei with StarDist", value=false) use_stardist
#@Float   (label="ROI specificity", value=2.0) roi_logpfa
#@Integer (label="Object channel", value=2) obj_channel
#@Float   (label="Object specificity", value=2.0) obj_logpfa
#@Float   (label="Downsampling", value=2,min=1,max=10,style="slider") downsampling_factor
#@Boolean (label="Remove ROI touching image borders", value=true) remove_border
#@Float   (label="Distance min [um]", value=0) distance_min
#@Float   (label="Distance max [um]", value=30) distance_max

/*
 * Analysis of the spatial distribution (spread) of markers
 * 
 * Objects in one channel are belonging to ROI defined in another channel.
 * 
 * If an image is opened: process the image
 * 
 * Input file: if a single image file process the image, if a csv file load condition and filename from the csv file
 * for batch processing. The csv file should have the columns condition, filename.
 * 
 * - Measure intensity for each ROI
 * - Measure spread of objects for each ROI
 *
 * Jerome Boulanger for Simon Bullock, Vanesa, Guilia (dec 2019 - 2021)
 */


var OUTOFBOUND = -100000; // flag value for out of bound
tbl1 = "Objects.csv";
tbl2 = "Regions of interest.csv";

starttime = getTime(); 

if (nImages !=0 ) {
	// First mode: process the opened image
	print("\n_____________________\nProcessing opened image");
	filename = getTitle();
	setBatchMode("hide");
	process_image(filename, condition, roi_channel, obj_channel, downsampling_factor);	
	setBatchMode("exit and display");
} else {
	if (filename.length() == 0) {
		print("\n_____________________\nProcessing a test image");
		generateTestImage();
		filename = getTitle();
		setBatchMode("hide");
		process_image(filename, condition, roi_channel, obj_channel, downsampling_factor);	
	setBatchMode("exit and display");
	} else if (endsWith(filename, ".csv")) {
		// Second mode: process the list of files
		print("\n_____________________\nProcessing list of files");		
		process_list_of_files(filename, roi_channel, obj_channel, downsampling_factor);		
	} else {
		// Third mode: process a single file
		print("\n_____________________\nProcessing single file");
		process_single_file(filename, condition, roi_channel, obj_channel, downsampling_factor);	
	}	
}

print("Processing finised in " + (getTime() - starttime)/60000 + "min");	

function process_list_of_files(tablename, roi_channel, obj_channel, downsampling_factor) {
	setBatchMode("hide");
	Table.open(tablename);
	N = Table.size;
	name = File.getName(tablename);
	for (n = 0; n < N; n++) {
		print("Image "+n+"/"+N);
		selectWindow(name);
		condition = Table.getString("condition", n);
		filename = Table.getString("filename", n);			
		print("Loading " + filename);
		run("Bio-Formats Importer", "open=["+filename+"] color_mode=Composite rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");	
		process_image(filename, condition, roi_channel, obj_channel, downsampling_factor);	
		close();
	}	
}

function process_single_file(filename, condition, roi_channel, obj_channel, downsampling_factor) {
	/* 
	 *  Process a single file and close the images at exit
	 */	
	print("Loading " + filename);
	run("Bio-Formats Importer", "open=["+filename+"] color_mode=Composite rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");	
	setBatchMode("hide");
	process_image(filename, condition, roi_channel, obj_channel, downsampling_factor);	
	setBatchMode("exit and display");
}

function process_image(filename, condition, roi_channel, obj_channel, downsampling_factor) {	
	closeRoiManager();	
	run("Stack to Hyperstack...", "order=xyczt(default) channels="+nSlices+" slices=1 frames=1 display=Color");	
	if (downsampling_factor != 1) {		
		w = round(getWidth()/downsampling_factor);
		h = round(getWidth()/downsampling_factor);
		run("Size...", "width="+w+" height="+h+" constrain average interpolation=None");		
	}
	id = getImageID;	
	if (use_stardist) {
		rois = segmentStarDist(getImageID, roi_channel);
	} else {
		rois = segment(getImageID, roi_channel, false, 20, roi_logpfa, 10, true);
	}
	objects = segment(getImageID, obj_channel, true, 5, obj_logpfa, 1, false);	
	dist = createSignedDistanceStack(id,rois);
	if (remove_border) {
		nrois = keepRoiNotTouchingImageEdges(rois,dist,10);
	} else {
		nrois = rois;
	}		
	//makeBandROI(dist, nrois, distance_min, distance_max);
	
	measureROIStats(id, dist, nrois, objects, obj_channel, condition, filename);

	
	selectImage(dist); close();
	selectImage(id);
	run("Make Composite");	
	roiManager("show none");
	Overlay.remove;
	addOverlays(nrois, true);
	addOverlays(objects, false);	
	closeRoiManager();		
}

function addOverlays(rois, addlabel) {	
	/*
	 * Add overlay for the roi whose indices are listed in the array 'rois'
	 */
	getPixelSize(unit, pixelWidth, pixelHeight);
	for (i = 0; i < rois.length; i++) {	
		roiManager("select", rois[i]);
		Overlay.addSelection(Roi.getStrokeColor, 1);
		Overlay.add();
		if (addlabel) {
			setColor("white");
			Overlay.drawString(""+rois[i]+1, getValue("X")/pixelWidth-5, getValue("Y")/pixelHeight);			
			Overlay.add();
		}
		Overlay.show();
	}
}

function generateTestImage() {
	newImage("Test image", "8-bit composite-mode", 500, 500, 3, 1, 1);
	run("Add Noise");
	run("Gaussian Blur...", "sigma=20");
	setSlice(1); run("Blue");
	makeOval(98, 84, 87, 85); run("Fill", "slice");	
	makeOval(264, 81, 130, 102); run("Fill", "slice");
	makeOval(118, 252, 119, 98); run("Fill", "slice");
	makeOval(321, 296, 129, 131); run("Fill", "slice");
	setSlice(2); run("Green");
	makeOval(407, 113, 14, 13); run("Fill", "slice");
	makeOval(406, 138, 7, 15); run("Fill", "slice");
	makeOval(388, 165, 13, 14); run("Fill", "slice");
	makeOval(366, 183, 12, 11); run("Fill", "slice");
	makeOval(390, 188, 6, 5); run("Fill", "slice");
	makeOval(418, 173, 7, 8); run("Fill", "slice");
	makeOval(423, 150, 10, 11); run("Fill", "slice");
	makeOval(202, 121, 6, 5); run("Fill", "slice");
	makeOval(197, 144, 4, 3); run("Fill", "slice");
	makeOval(182, 161, 9, 9); run("Fill", "slice");
	makeOval(166, 174, 2, 11); run("Fill", "slice");
	makeOval(137, 181, 10, 9); run("Fill", "slice");
	makeOval(189, 99, 11, 9); run("Fill", "slice");
	makeOval(253, 300, 20, 28); run("Fill", "slice");	
	makeOval(407, 369, 23, 24); run("Fill", "slice");
	makeOval(173, 52, 16, 4); run("Fill", "slice");
	makeOval(399, 70, 7, 17); run("Fill", "slice");
	makeOval(450, 112, 19, 16); run("Fill", "slice");
	makeOval(441, 166, 11, 9); run("Fill", "slice");
	makeOval(257, 96, 14, 9); run("Fill", "slice");
	makeOval(213, 150, 13, 7); run("Fill", "slice");
	makeOval(202, 174, 11, 9); run("Fill", "slice");
	makeOval(214, 131, 11, 7); run("Fill", "slice");
	makeOval(178, 185, 10, 6); run("Fill", "slice");
	run("Select None");
	run("Gaussian Blur...", "sigma=2");	
	run("Add Specified Noise...", "standard=2");
}

///////////////////////////////////////////////////////////////////////////////////
//                       Measurements
///////////////////////////////////////////////////////////////////////////////////
function createTable(name) {
	// creaate a table with title 'name' if it does not exist
	// return the size of the table
	if (!isOpen(name)) { Table.create(name);}
	selectWindow(name);
	return Table.size;
}

function makeBandROI(dist,rois, min, max) {
	// create a set of new ROI around each rois
	bands = newArray(rois.length);
	for (k = 0; k < rois.length; k++) {
		index = rois[k];
		roiManager("select", index);
		color = Roi.getStrokeColor();
		selectImage(dist);		
		setSlice(k+1);
		setThreshold(min, max);
		run("Create Selection");
		Roi.setStrokeColor(color);
		roiManager("update");		
		setSlice(k+1);
		run("Make Inverse");
		setColor(OUTOFBOUND);
		fill();
	}	
}

function measureRadialIntensities(id, dist, slices_index, dmin, dmax, n) {
	/* 
	 * Measure Mean Intensities in concentrics bands around ROIs
	 * Returns the mean for each bands in an array
	 * 
	 * TODO: take into account the objects masks
	 */
	vals = newArray(n);		
	for (i = 0; i < n; i++) {
		d0 = dmin + i / n * (dmax-dmin);
		d1 = dmin + i / n * (dmax-dmin);
		selectImage(dist);
		setSlice(slice_index);
		setThreshold(d0, d1);
		run("Create Selection");
		selectImage(id);
		run("Restore Selection");
		vals[i] = getValue("Mean");
	}
	return vals;
}

function getValues(id, slice, xpoints, ypoints) {
	vals = newArray(xpoints.length);
	selectImage(id);
	setSlice(slice);
	for (q = 0; q < xpoints.length; q++) {	
		vals[q] = getPixel(xpoints[q], ypoints[q]);
	}
	return vals;
}

function measureROIStats(id, dist, rois, objects, obj_channel, condition_name, image_name) {	
	N = 10;
	// compute object spread and other statistics for each ROI
	selectImage(id);
	getPixelSize(unit, pixelWidth, pixelHeight);
	// initialize accumulators
	swx = newArray(rois.length);
	swx2 = newArray(rois.length);
	swy = newArray(rois.length);
	swy2 = newArray(rois.length);	
	sw = newArray(rois.length);	 // sum of intensity in object per ROI
	count =  newArray(rois.length); // number of objects per ROI
	area =  newArray(rois.length); // area of the objects in the ROI	
	// allocate final measurements arrays
	spread = newArray(rois.length); 
	entropy = newArray(rois.length); // sum[ w/(sum w) log(w/(sum w))] = (sum(w log w) - (sum w) log(sum w)) / sum w
	distance = newArray(rois.length); // avergae object distance to the  nuclei
	radial_intensity =  newArray(N * rois.length); // radial intensities
	
	for (k = 1; k <= nSlices; k++) { // index on the distance stack						
		for (i = 0; i < objects.length; i++) {
			selectImage(dist);
			setSlice(k);
			roiManager("select", objects[i]);
			val = getValue("Min");			
			// if the minimum distance in the object is more than OUTOFBOUND, the object is considered for analysis
			if (val > OUTOFBOUND) {
				count[k-1] = count[k-1] + 1;				
				distance[k-1] = distance[k-1] + getValue("Mean");				
				// get the coordinates of points inside the object ROI
				roiManager("select", objects[i]);				
				Roi.getContainedPoints(xpoints, ypoints);
				// get the distances values
				di = getValues(dist, k, xpoints, ypoints);
				
				// get the intensity values
				w = getValues(id, obj_channel, xpoints, ypoints);				
				// compute all accumulators
				for (q = 0; q < xpoints.length; q++) {									
					dx = xpoints[q] * pixelWidth;
					dy = ypoints[q] * pixelHeight;
					sw[k-1]      = sw[k-1]   + w[q];
					swx[k-1]     = swx[k-1]  + w[q] * dx;
					swy[k-1]     = swy[k-1]  + w[q] * dy;
					swx2[k-1]    = swx2[k-1] + w[q] * dx * dx;
					swy2[k-1]    = swy2[k-1] + w[q] * dy * dy;
					if (w[q] > 0.1) {
						entropy[k-1] = entropy[k-1] + w[q] * log(w[q]);
					}
					area[k-1]    = area[k-1] + 1;				
					// compute the index for radial intensity array										
					dindex = round((di[q] * pixelWidth - distance_min) / (distance_max -distance_min) * (N-1));
					if (dindex >= 0 && dindex < N) {
						radial_intensity[dindex + N * (k-1)] = radial_intensity[dindex + N * (k-1)] + w[q];
					}
				}				
			}
		}
		// derive spread, entropy, distances, etc from the accumulators
		swx[k-1] = swx[k-1] / sw[k-1];
		swy[k-1] = swy[k-1] / sw[k-1];
		swx2[k-1] = swx2[k-1] / sw[k-1] - swx[k-1] * swx[k-1];
		swy2[k-1] = swy2[k-1] / sw[k-1] - swy[k-1] * swy[k-1];
		spread[k-1] = sqrt(swx2[k-1] + swy2[k-1]);	
		entropy[k-1] = (sw[k-1]*log(sw[k-1]) - entropy[k-1]) / sw[k-1];
		distance[k-1] = distance[k-1] / count[k-1] * pixelWidth;
		// normalize the radial intensity by the area of the band
		selectImage(dist);
		setSlice(k);
		for (dindex = 0; dindex < N; dindex++) {
			lower = (distance_min + dindex * (distance_max - distance_min) / (N-1)) / pixelWidth;
			upper = (distance_min + (dindex+1) * (distance_max - distance_min) / (N-1)) / pixelWidth;
			setThreshold(lower, upper);
			run("Create Selection");
			a = getValue("Area");
			radial_intensity[dindex + N * (k-1)] = radial_intensity[dindex + N * (k-1)] / a;
		}
	}
	selectImage(id);
	Stack.getDimensions(width, height, channels, slices, frames);
	roitbl = "Regions of interest";
	n = createTable(roitbl);
	channel_names = split(channel_names_str,',');	
	for (k = 0; k < rois.length; k++) {
		if (count[k] > 0) {					
			index = rois[k];			
			roiManager("select", index);
			Table.set("Condition"    , n, condition_name);
			Table.set("Image"        , n, File.getName(image_name));
			Table.set("ROI channel"  , n, channel_names[roi_channel-1]);
			Table.set("Object channel"  , n, channel_names[obj_channel-1]);
			Table.set("Image"        , n, File.getName(image_name));
			Table.set("roiID"        , n, index+1);	
			for (c = 1; c <= channels; c++) {
				Stack.setChannel(c);
				Table.set("Mean ["+channel_names[c-1]+"]", n, getValue("Mean"));
			}
			Table.set("Object Area"  , n, area[k]);
			Table.set("Object Count"  , n, count[k]);
			Table.set("Object Spread", n, spread[k]);
			Table.set("Object Entropy", n, entropy[k]);
			Table.set("Object Mean Intensity", n, sw[k]/area[k]);
			Table.set("Object Mean Distance", n, distance[k]);	
			for (dindex = 0; dindex < N; dindex++) {
				dvalue = round(distance_min + dindex * (distance_max - distance_min) / (N-1));
				Table.set("Intensity at Distance " + dvalue +" "+ unit, n, radial_intensity[dindex + N * k]);
			}		
			n++;
		}
	}	
	Table.update();
}

///////////////////////////////////////////////////////////////////////////////////
//                       Distance functions
///////////////////////////////////////////////////////////////////////////////////
function createSignedDistanceStack(id,rois) {
	print(" - Compute distance map");
	getDimensions(width, height, channels, slices, frames);
	newImage("SignedDistanceStack", "32-bit black", width, height, rois.length);
	stack = getImageID();
	for (i = 0; i < rois.length; i++) {
		d = createSignedDistance(id,rois[i]);
		run("Copy");
		selectImage(stack);
		setSlice(i+1);		
		run("Paste");	
		selectImage(d); close();
	}
	selectImage(stack);
	run("Select None");
	for (y = 0; y < height; y++) {
		for (x = 0; x < width; x++) {	
			idx = -1;
			min = 1e6;
			for (n = 1; n <= rois.length; n++) {
				setSlice(n);		
				val = getPixel(x,y);
				if (val < min) {
					min = val;
					idx = n;
				}
			}
			for (n = 1; n <= rois.length; n++) {				
				if (n != idx ) {
					setSlice(n);		
					setPixel(x,y,OUTOFBOUND); 
				} 
			}
		}
	}	
	N = nSlices;
	if (N != rois.length) {
		print("\nERROR in createSignedDistanceStack");
		print("  >> The number of slices in the distance image and of ROI in the list are different.\n");
	}
	return stack;
}

function createSignedDistance(id,roi) {
	run("Options...", "iterations=1 count=1 black edm=32-bit");
	selectImage(id);
	getDimensions(width, height, channels, slices, frames);
	newImage("Untitled", "8-bit black", width, height, 1);
	in0 = getImageID;
	roiManager("select", roi);
	setColor(255, 255, 255);
	fill();	
	run("Distance Map");
	in = getImageID;
	run("Multiply...", "value=-1");
	newImage("Untitled", "8-bit white", width, height, 1);
	out0 = getImageID;
	roiManager("select", roi);
	setColor(0, 0, 0);
	fill();	
	run("Distance Map");
	out = getImageID;
	imageCalculator("Add create 32-bit", in, out);
	dist = getImageID;
	selectImage(in0); close();
	selectImage(in); close();
	selectImage(out0); close();
	selectImage(out); close();
	selectImage(dist);	
	return dist;
}


function keepRoiNotTouchingImageEdges(roi, dist, b) {
	/*
	 * Keep roi close to image edge and return the list of index
	 * 
	 * does not modify the ROI manager but remove the frames from distance image dist
	 * 
	 */
	 
	selectImage(dist);
	N = nSlices;
	if (N != roi.length) {
		print("\nERROR in keepRoiNotTouchingImageEdges");
		print(" >> The number of slices in the distance image and of ROI in the list are different.\n");
	}	
	flag = newArray(roi.length);// flag true if the ROI is not touching the edge
	w = getWidth();
	h = getHeight();
	makeRectangle(b, b, w-2*b, h-2*b);
	run("Make Inverse");
	Roi.getContainedPoints(xpoints, ypoints);
	K = 0; // count number of region inside
	for (n = 0;  n < roi.length; n++) {
		selectImage(dist);	
		setSlice(n+1);		
		flag[n] = true;		
		for (i = 0; i < xpoints.length; i++) {			
			if (getPixel(xpoints[i],ypoints[i]) <  b && getPixel(xpoints[i],ypoints[i]) != OUTOFBOUND)  {
				flag[n] = false;
			}			
		}
		if (flag[n] == true) {
			K = K + 1; // count the number of ROI to keep
		}
	}	
	// create the list of ROI to keep
	nroi = newArray(K); 
	k = 0;
	for (n = 0; n < roi.length; n++) {		
		if (flag[n] == true) {
			nroi[k] = roi[n];
			k = k + 1;
		} else {
			//print("Removing " + (n+1));
			roiManager("select", roi[n]);
			Roi.setStrokeColor("#FF0000");
		}
	}
	// remove image from the distance stack
	for (n = N-1;  n >= 0; n--) {
		if (flag[n] == false) {
			setSlice(n+1);
			if (nSlices > 1) {
				run("Delete Slice");			
			} else {
				run("Select All");
				run("Cut");
				setColor(OUTOFBOUND, OUTOFBOUND, OUTOFBOUND);
				fill();
			}
		}
	}
	print("- Keeping " + nroi.length + " / " + roi.length + " objects away from image border");
	selectImage(dist);	
	if (nSlices != nroi.length) {
		print("\nERROR in keepRoiNotTouchingImageEdges (end)");
		print("    The number of slices in the distance image and of ROI in the list are different.");
		print("    nroi.length : " +  nroi.length );
		print("    nSlices     : " +  nSlices);
	}
	return nroi;
}

///////////////////////////////////////////////////////////////////////////////////
//                     Segmentation function
///////////////////////////////////////////////////////////////////////////////////

function segmentStarDist(id, ch) {
	// Segment the image with StarDist and returns the list of ROIs
	print(" - segment  channel #" + ch + " with StarDist");
	selectImage(id);
	run("Select None");		
	run("Duplicate...", "title=[mask-"+ch+"] duplicate channels="+ch);
	id1 = getImageID();
	n0 = roiManager("count");	
	run("Command From Macro", "command=[de.csbdresden.stardist.StarDist2D], args=['input':'mask-"+ch+"', 'modelChoice':'Versatile (fluorescent nuclei)', 'normalizeInput':'true', 'percentileBottom':'1.0', 'percentileTop':'99.8', 'probThresh':'0.5', 'nmsThresh':'0.4', 'outputType':'ROI Manager', 'nTiles':'1', 'excludeBoundary':'2', 'roiPosition':'Automatic', 'verbose':'false', 'showCsbdeepProgress':'false', 'showProbAndDist':'false'], process=[false]");
	n1 = roiManager("count");
	a = generateSequence(n0,n1-1);
	cmap = newArray("#5555AA","#55AA55","#FF9999","#FFFF55");
	for (i = n0; i < n1; i++) {
		roiManager("select", i);
		Roi.setStrokeColor(cmap[ch-1]);
	}
	selectImage(id1); close();
	selectImage(id);
	roiManager("show all without labels");
	print("   - "+(a.length)+" ROI added");
	return a;
}

function segment(id, ch, sharpen, bg, logpfa, minsz, fholes) {
	/* segment the channel ch on image image id	
	 *  id: image ID
	 *  ch : channel
	 *  sharpen: sharpen the image first
	 *  bg : size of the rolling ball in um
	 *  logpfa : - log10 probability of false alarm
	 *  minsz : typical diameter
	 *  fhole : fill holes
	 */
	getPixelSize(unit, px, py);
	r0 = round(1/px); 
	r1 = 0.25 / px;
	print(" - segment channel #" + ch + " ( bg = "+bg+", minsize = "+minsz+", fill = "+fholes+" )");
	selectImage(id);
	run("Select None");		
	run("Duplicate...", "title=[mask-"+ch+"] duplicate channels="+ch);	
	run("32-bit");
	id1 = getImageID();
	getPixelSize(unit, pw, ph);
	if (sharpen) {		
		run("Gaussian Blur...", "sigma="+r1);
		run("Convolve...", "text1=[-0.16 -0.66 -0.16\n -0.66 5 -0.66\n-0.16 -0.66 -0.16\n] ");			
	} 
	run("Median...", "radius="+r0);	
	run("Subtract Background...", "rolling="+(bg/px)+" disable");
	setImageThreshold(logpfa);
	run("Convert to Mask");
	if (fholes) {
		run("Fill Holes");
		for (i = 0; i < 5; i++) {
			run("Maximum...", "radius="+r0);
			run("Minimum...", "radius="+r0);
			run("Median...", "radius="+r0);
		}
	}
	run("Watershed");
	n0 = roiManager("count");
	run("Analyze Particles...", "size="+(PI* minsz*minsz/4)+"-Infinity add");
	n1 = roiManager("count");
	a = generateSequence(n0,n1-1);
	// assign a color to the ROI according to the channel
	cmap = newArray("#5555AA","#55AA55","#FF9999","#FFFF55");
	for (i = n0; i < n1; i++) {
		roiManager("select", i);
		Roi.setStrokeColor(cmap[ch-1]);
	}
	selectImage(id1); close();
	print("   - "+(n1-n0+1)+" ROI added");
	selectImage(id);
	roiManager("show all without labels");
	return a;
}

function generateSequence(a,b) {	
	// create an array with a sequence from a to b with step 1
	A = newArray(b-a+1);
	for (i = 0; i < A.length; i++) {A[i] = a+i;}
	return A;
}

function setImageThreshold(logpfa) {	
	// Compute a threshold from a quantile on robust statistics
	A = convertImageToArray();
	med = computeArrayQuantile(A,0.1);
	mad = computeMeanAbsoluteDeviation(A,med);		
	t = med + logpfa * log(10) * mad;
	print("  - med = " + med  + " mad = " + mad + " threshold = " + t);
	Array.getStatistics(A, min, max);
	setThreshold(t, max);
	return t;
}

function computeArrayQuantile(A,q) {
	P = Array.rankPositions(A);	
	return A[P[round(q*A.length)]];
}

function computeArrayMedian(A) {
	P = Array.rankPositions(A);	
	return A[P[round(A.length/2.0)]];
}

function computeMeanAbsoluteDeviation(A, med) {
	s = 0;
	for (i = 0; i < A.length; i++) {
		s += abs(A[i] - med);
	}
	return s / A.length;
}

function convertImageToArray() {
	/* 
	 *  Convert the current image to an array
	 */
	width = getWidth();
	height = getHeight();
	A = newArray(width*height);	
	for (y = 0; y < height; y++) {
		for (x = 0; x < width; x++) {
			A[x + width * y] = getPixel(x, y);
		}
	}
	return A;
}

///////////////////////////////////////////////////////////////////////////////////
// Utilities
///////////////////////////////////////////////////////////////////////////////////

function closeRoiManager() {
	closeWindow("ROI Manager");
}

function closeWindow(name) {
	if (isOpen(name)) {
		selectWindow(name);
		run("Close");
	}
}


