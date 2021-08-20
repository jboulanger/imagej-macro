#@File    (label="Input File",description="Either a csv file or an image file") filename
#@String  (label="Group",description="Experimental condition", value="control") condition
#@Integer (label="ROI channel", value=1) roi_channel
#@Integer (label="Object channel", value=2) obj_channel
#@String  (label="Channels name", value="DAPI,GFP") channel_names_str
#@Float   (label="Downsampling", value=2,min=1,max=10,style="slider") downsampling_factor
#@Boolean (label="Remove ROI touching image borders", value=true) remove_border
#@Float   (label="Distance min", value=0) distance_min
#@Float   (label="Distance max", value=100) distance_max

/*
 * Analysis of the spatial distribution (spread) of markers
 * 
 * Works for multi channel images
 * 
 * Input file: if a single image file process the image, if a csv file load condition and filename from the csv file
 * for batch processing.
 * 
 * Measure intensity for each ROI
 *
 * Jerome Boulanger for Simon Bullock, Vanesa, Guilia (dec 2019 - 2021)
 */


var OUTOFBOUND = -100000; // flag value for out of bound
tbl1 = "Objects.csv";
tbl2 = "Regions of interest.csv";

starttime = getTime(); 

if (nImages !=0 ) {
	// First mode: process the opened image
	print("Processing opened image");
	filename = getTitle();
	process_image(filename, condition, roi_channel, obj_channel, downsampling_factor);
} else {
	if (endsWith(filename, ".csv")) {
		// Second mode: process the list of files
		print("Processing list of files");
		setBatchMode("hide");
		process_list_of_files(filename, roi_channel, obj_channel, downsampling_factor);
		setBatchMode("exit and display");
	} else {
		// Third mode: process a single file
		print("Processing single file");
		process_single_file(filename, condition, roi_channel, obj_channel, downsampling_factor);
	}	
}

print("Processing finised in " + (getTime() - starttime)/60000 + "min");	


function process_list_of_files(tablename, roi_channel, obj_channel, downsampling_factor) {
	Table.open(tablename);
	N = Table.size;
	name = File.getName(tablename);
	for (n = 0; n < N; n++) {
		print("Image "+n+"/"+N);
		selectWindow(name);
		condition = Table.getString("condition", n);
		filename = Table.getString("filename", n);	
		process_single_file(filename, condition, roi_channel, obj_channel, downsampling_factor);
	}
}

function process_single_file(filename, condition, roi_channel, obj_channel, downsampling_factor) {
	/* Process a single file and close the images at exit
	 * 
	 */	
	print("Loading " + filename);
	run("Bio-Formats Importer", "open=["+filename+"] color_mode=Composite rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
	process_image(filename, condition, roi_channel, obj_channel, downsampling_factor);
	run("Close All");
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

	rois =    segment(id, roi_channel, false, 100, 0.01, 1000,true);
	objects = segment(id, obj_channel, true , 10 , 0.00001, 1,false);
	setBatchMode(true);
	dist = createSignedDistanceStack(id,rois);
	if (remove_border) {
		nrois = keepRoiNotTouchingImageEdges(rois,dist,10);
	} else {
		nrois = rois;
	}	
	setBatchMode(true);
	makeBandROI(dist, nrois, distance_min, distance_max);
	measureROIStats(id, dist, rois, objects, obj_channel, condition, filename);
	selectImage(id);
	run("Make Composite");
	roiManager("show all without labels");	
}

///////////////////////////////////////////////////////////////////////////////////
//                       Measurements
///////////////////////////////////////////////////////////////////////////////////
function createTable(name) {
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
		selectImage(dist);		
		setSlice(k+1);
		setThreshold(min, max);
		run("Create Selection");
		roiManager("update");
		setSlice(k+1);
		run("Make Inverse");
		setColor(OUTOFBOUND);
		fill();
	}	
}

function measureROIStats(id, dist, rois, objects, obj_channel, condition_name, image_name) {
	
	// compute spread and other statistics on objects mask
	swx = newArray(rois.length);
	swx2 = newArray(rois.length);
	swy = newArray(rois.length);
	swy2 = newArray(rois.length);
	sw = newArray(rois.length);	 // sum of intensity in object per ROI
	count =  newArray(rois.length); // number of objects per ROI
	area =  newArray(rois.length); // area of the objects in the ROI	
	spread = newArray(rois.length); 
	entropy = newArray(rois.length); // sum[ w/(sum w) log(w/(sum w))] = (sum(w log w) - (sum w) log(sum w)) / sum w
	
	
	
	for (k = 1; k <= nSlices; k++) { // index on the distance stack						
		for (i = 0; i < objects.length; i++) {
			selectImage(dist);
			setSlice(k);	
			roiManager("select", objects[i]);
			val = getValue("Mean");
			if (abs(val - OUTOFBOUND) > 1) {
				count[k-1] = count[k-1] + 1;
				selectImage(id);
				roiManager("select", objects[i]);
				Stack.setChannel(obj_channel);
				Roi.getContainedPoints(xpoints, ypoints);
				for (q = 0; q < xpoints.length; q++) {					
					w = getPixel(xpoints[q], ypoints[q]);
					dx = xpoints[q];
					dy = ypoints[q];
					sw[k-1]      = sw[k-1]   + w;
					swx[k-1]     = swx[k-1]  + w * dx;
					swy[k-1]     = swy[k-1]  + w * dy;
					swx2[k-1]    = swx2[k-1] + w * dx * dx;
					swy2[k-1]    = swy2[k-1] + w * dy * dy;
					if (w > 1) {
						entropy[k-1] = entropy[k-1] + w * log(w);
					}
					area[k-1]    = area[k-1] + 1;
				}
			}
		}
		swx[k-1] = swx[k-1] / sw[k-1];
		swy[k-1] = swy[k-1] / sw[k-1];
		swx2[k-1] = swx2[k-1] / sw[k-1] - swx[k-1] * swx[k-1];
		swy2[k-1] = swy2[k-1] / sw[k-1] - swy[k-1] * swy[k-1];
		spread[k-1] = sqrt(swx2[k-1] + swy2[k-1]);	
		entropy[k-1] = (sw[k-1]*log(sw[k-1]) - entropy[k-1]) / sw[k-1];
	}
	selectImage(id);
	Stack.getDimensions(width, height, channels, slices, frames);
	roitbl = "Regions of interest";
	n = createTable(roitbl);
	for (k = 0; k < rois.length; k++) {		
		index = rois[k];
		roiManager("select", index);
		Table.set("Condition"    , n, condition_name);
		Table.set("Image"        , n, File.getName(image_name));
		Table.set("roiID"        , n, index+1);	
		for (c = 1; c <= channels; c++) {
			Stack.setChannel(c);
			Table.set("Mean Ch"+c, n, getValue("Mean"));
		}
		Table.set("Object Area"  , n, area[k]);
		Table.set("Object Count"  , n, count[k]);
		Table.set("Object Spread", n, spread[k]);
		Table.set("Object Entropy", n, entropy[k]);
		Table.set("Object Mean Intensity", n, sw[k]/area[k]);
		n++;
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
function segment(id, ch, sharpen, bg, pfa, minsz, fholes) {
	// segment the channel ch on image image id
	print(" - segment channel #" + ch + " ( bg = "+bg+", minsize = "+minsz+", fill = "+fholes+" )");
	selectImage(id);
	run("Select None");		
	run("Duplicate...", "title=[mask-"+ch+"] duplicate channels="+ch);	
	run("32-bit");
	id1 = getImageID();
	if (sharpen) {		
		run("Square Root");
		run("Gaussian Blur...", "sigma=1");
		run("Convolve...", "text1=[-0.16 -0.66 -0.16\n -0.66 4 -0.66\n-0.16 -0.66 -0.16\n] ");
		
	} else {
		run("Median...", "radius=3");	
	}	
	run("Subtract Background...", "rolling="+bg+" disable");	
	getStatistics(area, mean, min, max, std);
	t = computeImageThreshold(pfa);	
	print("   - threshold  with pfa " + pfa + " : " + t);
	setThreshold(t, max);		
	run("Convert to Mask");
	if (fholes) {
		run("Open");
		run("Fill Holes");
		for (i = 0; i < 3; i++) {	
			//run("Open");
			run("Median...", "radius=5");
		}
	}
	run("Watershed");	
	n0 = roiManager("count");	
	run("Analyze Particles...", "size="+minsz+"-Infinity pixel add");		
	n1 = roiManager("count");
	a = generateSequence(n0,n1-1);
	// assign a color to the ROI according to the channel
	cmap = newArray("#5555AA","#55AA55","#FF9999","#FFFF55");
	for (i = n0; i < n1; i++) {	
		roiManager("select", i)	;
		Roi.setStrokeColor(cmap[ch-1]);
	}
	selectImage(id1); close();
	print("   - "+(n1-n0+1)+" ROI added");
	return a;	
}

function generateSequence(a,b) {	
	// create an array with a sequence from a to b with step 1
	A = newArray(b-a+1);
	for (i = 0; i < A.length; i++) {A[i] = a+i;}
	return A;
}

function computeImageThreshold(pfa) {	
	// Compute a threshold from a quantile on robust statistics
	getStatistics(area, mean, min, max, std, histogram);
	med = 0;
	n = 0;
	while (n < area/2) {
		n = n + histogram[med];
		med = med + 1;
	}
	med = med * (max - min) / histogram.length + min;
	run("Duplicate...", "title=tmp ");
	tmp = getImageID();
	run("32-bit");
	run("Subtract...", "value="+med);
	run("Abs");
	getStatistics(area, mean, min, max, std, histogram);
	mad = 0;	
	n = 0;
	while (n < area/2) {
		n = n + histogram[med];
		mad = mad + 1;
	}
	mad = mad * (max - min) / histogram.length + min;
	mad = 1.48 * mad;		
	selectImage(tmp); close();
	print("med = " + med  + " mad = " + mad);
	return med - 2 * norminv(pfa) * mad;
}

function norminv(x) {
	if (x < 0.25) {
		return -sqrt(  log(1/(x*x)) - 0.9*log(log(1/(x*x))) - log(2*3.14159) );
	} else { 
		return 0;
	}
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


