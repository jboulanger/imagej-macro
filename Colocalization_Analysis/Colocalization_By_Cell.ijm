#@File    (label="Input file",description="Either a csv file or an image file, tye test to generate image") filename
#@String  (label="Experimental group",description="Experimental group", value="control") condition
#@String  (label="Channels name", description="List the name of the channels in order separated by commas", value="DAPI,GFP") channel_names_str
#@Integer (label="ROI channel", description="Indicate the index of the channel used to segment ROI (eg nuclei)", value=1) roi_channel
#@String  (label="Obj channels", description="List the indices of the channels used for colocalization", value="2,3") obj_channels
#@String  (label="Mask channel", description="Indicate if necessary the index of the channel used for masking cells", value=2) mask_channel
#@Boolean (label="Mask background", description="Tick this to limit the ROI by the signal in the mask channel", value=true) maskbg
#@Float   (label="Downsampling", description="Downsampling factor", value=2,min=1,max=10,style="slider") downsampling_factor
#@Boolean (label="Segment ROI with StarDist", description="Indicate whether you want to use startdist to segment nuclei", value=false) use_stardist
#@Float   (label="ROI specificity", description="If not using stardist, define a speficity level to define the threshold", value=2.0) roi_logpfa
#@Boolean (label="Remove ROI touching image borders", description="Tick this to remove cell touching image border", value=true) remove_border
#@String  (label="preprocessing", description="type of object detection", choices={"LoG","Rolling Ball","None"}) preprocess_type
#@Float   (label="scale [px]", description="scale of the object to detect", value=2, min=1, max=10, style="slider") scale
#@Float   (label="Obj Specificity",description="set the specificity to segmenting objects", value=0.25) obj_logpfa
#@Boolean (label="Colocalize on original?", description="use the original image of the preprocessed image", value=true) onoriginal
#@String  (label="overlay colors", description="comma separated list of colors for control overlay", value="blue,red,green,white") colors_str
#@String  (label="Table name", value="colocalization.csv") tbl1

/*
 * Analysis of the colocalization of markers cell by cell
 * 
 * Objects in one channel are belonging to ROI defined in another channel.
 * 
 * If an image is opened: the macro will process the image
 * 
 * Input file: if a single image file process the image, if a csv file load condition and filename from the csv file
 * for batch processing. The csv file should have the columns condition, filename.
 * 
 * - Measure intensity for each ROI
 * - Measure spread of objects for each ROI
 *
 * Jerome Boulanger for Vanesa (2021)
 */

var OUTOFBOUND = -100000; // flag value for out of bound

var starttime = getTime(); 

print("\\Clear");

// convert the comma-separted string of channels to an array
channel_names = split(channel_names_str,',');
obj_channels = parseCSVnumbers(obj_channels);

// inialize default keys, vals list to empty array
keys = newArray("ROI Specificity", "Object Specificity", "Zoom", "Remove Border", "Mask Background");
vals = newArray(roi_logpfa, obj_logpfa, downsampling_factor, remove_border, maskbg);

if (nImages !=0 ) {
	// First mode: process the opened image
	print("\n_____________________\nProcessing opened image");
	filename = getTitle();	
	setBatchMode("hide");
	process_image(filename, condition, channel_names, roi_channel, obj_channels, mask_channel, keys, vals);	
	setBatchMode("exit and display");
} else {
	if (matches(File.getName(filename),"test")) {
		print("\n_____________________\nProcessing a test image");
		generateTestImage();
		filename = getTitle();
		setBatchMode("hide");
		roi_channel = 1;		
		downsampling_factor = 1;
		channel_names = newArray("Blue","Red","Green");
		obj_channels = newArray(2,3);
		mask_channel = 0;		
		process_image(filename, condition, channel_names, roi_channel, obj_channels, mask_channel, keys, vals);
		setBatchMode("exit and display");
	} else if (endsWith(filename, ".csv")) {
		// Second mode: process the list of files
		print("\n_____________________\nProcessing list of files " + filename);		
		process_list_of_files(filename, condition, channel_names, roi_channel, obj_channels, mask_channel, keys,vals);		
	} else {
		// Third mode: process a single file
		print("\n_____________________\nProcessing single file");		
		process_single_file(filename, condition, channel_names, roi_channel, obj_channels, mask_channel, keys, vals);	
	}
}

print("Processing finised in " + (getTime() - starttime)/60000 + "min");	


function parseCSVnumbers(str) {
	// return array of numbers from comma separated string
	a = split(str,",");
	b = newArray(a.length);
	for (i=0;i<a.length;i++){b[i] = parseFloat(a[i]);}
	return b;
}


function getColors(str) {
	colors = split(str,",");
	Stack.getDimensions(width, height, channels, slices, frames);
	while (colors.length < channels) {
		colors = Array.concat(colors, "yellow");
	}
	return colors;
}

function process_list_of_files(tablename, condition, channel_names, roi_channel, obj_channel, mask_channel, keys, vals) {
	/*
	 * Process a files listed in a csv file 'tablename'
	 * The table can have the columns 'Filename', 'Condition', 'Channel Names', 'ROI channel', 'Object Channel'
	 * 
	 * Channel Names are coma separated names between double quotes : "A,B,C,D"
	 */
	setBatchMode("hide");
	Table.open(tablename);	
	N = Table.size;
	name = File.getName(tablename);
	path = File.getDirectory(tablename);
	hasCondition = isColumnInTable("Condition");
	hasChannelNames = isColumnInTable("Channel Names");
	hasRoi = isColumnInTable("ROI Channel");
	hasObj = isColumnInTable("Object Channel");
	hasMask = isColumnInTable("Mask Channel");
	
	new_keys = listOtherColumns(newArray("Filename","Condition","Channel Names","ROI Channel","Object Channel","Mask Channel"));
	new_vals = newArray(new_keys.length);
	
	for (n = 0; n < N; n++) {
		if (n==0) {
			print("Image "+(n+1)+"/"+N + "");
		} else {
			rtime = (N-n) * (getTime() - starttime) / n;
			print("Image "+(n+1)+"/"+N + " remaining time " + (rtime / 60000)+"min");
		}
		selectWindow(name);
		filename = Table.getString("Filename", n);
		if (hasCondition) {
			condition = Table.getString("Condition", n);		
		}
		if (hasChannelNames) {
			channel_names_str = Table.getString("Channel Names", n);			
			channel_names = split(channel_names_str,",");			
		}
		if (hasRoi) {
			roi_channel = Table.getString("ROI Channel", n);
		}
		if (hasObj) {
			obj_channels = split(Table.getString("Object Channel", n), ",");
		}
		if (hasMask) {
			mask_channel = Table.getString("Mask Channel", n);
		}
		for (i = 0; i < new_keys.length; i++) {
			new_vals[i] = Table.getString(new_keys[i], n);
		}
		all_keys = Array.concat(new_keys, keys);	
		all_vals = Array.concat(new_vals, vals);		
		
		print("Loading " + filename + " with condition " + condition);		
		run("Bio-Formats Importer", "open=["+path + File.separator + filename+"] color_mode=Composite rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");			
		process_image(filename, condition, channel_names, roi_channel, obj_channels, mask_channel, all_keys, all_vals);		
		close();
		print("\n");
		// save temporary version of the roi and obj table in case of crash.
		selectWindow(tbl1);
		Table.save(path + File.separator + "roi-tmp.csv");
		selectWindow(tbl2);
		Table.save(path + File.separator + "obj-tmp.csv");		
		selectWindow("Log");
		saveAs("Text", path + File.separator + "log.txt");		
	}	
}

function isColumnInTable(colname) {
	/* Return true if the table has a column with name 'colname' */
	headers = split(Table.headings(),'\t');
	for (k = 0; k < headers.length; k++) {
		if (matches(headers[k], colname)) {
			return true;
		}
	}
	return false;
}

function listOtherColumns(colnames) {
	/* Return an array with the columns name which are not in colnames	 */
	headers = split(Table.headings(),'\t');
	other = headers;
	j = 0;
	for (k = 0; k < headers.length; k++) {
		keep = true;
		for (i = 0; i < colnames.length && keep; i++) {
			if (matches(headers[k], colnames[i])) {
				keep = false;
			}
		}
		if (keep) {
			other[j] = headers[k];
			j = j + 1;
		}
	}	
	return Array.trim(other,j);
}

function process_single_file(filename, condition, channel_names, roi_channel, obj_channel, mask_channel, keys, vals) {
	/* 
	 *  Process a single file and close the images at exit
	 */	
	print("Loading " + filename);
	run("Bio-Formats Importer", "open=["+filename+"] color_mode=Composite rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");	
	setBatchMode("hide");
	process_image(filename, condition, channel_names, roi_channel, obj_channel, mask_channel, keys, vals);	
	setBatchMode("exit and display");
}

//function process_image(filename, condition, channel_names, roi_channel, obj_channel, mask_channel, downsampling_factor, keys, vals) {
function process_image(filename, condition, channel_names, roi_channel, obj_channels, mask_channel, keys, vals) {
	/*	
	 * 	Process the opened image
	 * 	- resample the image
	 * 	- segment ROI in roi_channel
	 * 	- compute distance map 
	 */	

	keys = Array.concat(keys, newArray("Filename","Condition"));
	vals = Array.concat(vals, newArray(filename,condition));
	
	colors = getColors(colors_str); 
	chstr = "Channel Names : ";
	for (i = 0; i < channel_names.length; i++) {
		channel_names[i] = String.trim(channel_names[i]);
		channel_names[i] = replace(channel_names[i],"\"","");
		chstr += "[" + toString(i+1) + ":" +channel_names[i] + "]";
		if (i != channel_names.length-1) {chstr+=", ";}
	}	
	print(chstr);	

	closeRoiManager();	
	run("Stack to Hyperstack...", "order=xyczt(default) channels="+nSlices+" slices=1 frames=1 display=Color");	
	if (downsampling_factor != 1) {
		print(" - Downsampling by a "+ downsampling_factor + "x");
		w = round(getWidth()/downsampling_factor);
		h = round(getWidth()/downsampling_factor);
		run("Size...", "width="+w+" height="+h+" constrain average interpolation=None");		
	}
	img1 = getImageID;
	
	Stack.getDimensions(width, height, channels, slices, frames);
	if (channel_names.length != channels) {
		exit("*** Error not matching channel names with number of channels !***");
		return 1;
	}
	
	if (use_stardist) {
		rois = segmentStarDist(getImageID, roi_channel);
	} else {
		rois = segment(getImageID, roi_channel, false, 50, roi_logpfa, 10, true);
	}
	
	dist = createSignedDistanceStack(img1, rois);	
	
	if (remove_border) {
		nrois = keepRoiNotTouchingImageEdges(rois,dist,10);
	} else {
		nrois = Array.copy(rois);
	}
	
	if (maskbg) {
		bg = segmentBackground(img1, mask_channel, downsampling_factor);
		applyBackground(dist, bg);	
		selectImage(bg); close();
	}
	closeRoiManager();	
	parents = getDistROI(dist);	
	selectImage(dist); close();				
	selectImage(img1);run("Select None");
	img2 = preprocess(img1, preprocess_type, scale);
	print("Number of Parent ROI: " + parents.length);
	thresholds = getChannelThresholds(img2, obj_channels, obj_logpfa);		
	children = getAllchildrenROI(img2, parents, obj_channels, thresholds, "0-Infinity", "0.0-1.0", colors);
	print("Number of Children ROI: " + children.length);	
	children_channel = getROIChannels(children);
	measureColocalization(tbl1, img1, img2, obj_channels, thresholds, parents, children, children_channel, onoriginal, channel_names, keys,vals);

	selectImage(img2);close();
	selectImage(img1);
	run("Make Composite");	
	roiManager("show none");
	Overlay.remove;	
	addOverlays(parents, true);
	addOverlays(children, false);		
	closeRoiManager();	
	run("Collect Garbage");
	run("Select None");	
	return 0;
}


function applyBackground(dist, bg) {
	/*
	 * Set distance stack values to outofbound if in background
	 */
	selectImage(bg);
	run("Create Selection");	
	run("Make Inverse");
	selectImage(dist);
	setColor(OUTOFBOUND);
	for (k = 1 ; k <= nSlices; k++) {
		setSlice(k);
		run("Restore Selection");
		fill();
	}
}

function removeROI(rois) {
	a = Array.copy(rois);
	Array.sort(a);
	for (i = a.length-1; i >= 0; i--) {
		roiManager("select", a[i]);
		roiManager("delete");
	}
}

function addOverlays(rois, addlabel) {	
	/*
	 * Add overlay for the roi whose indices are listed in the array 'rois'
	 */
	getPixelSize(unit, pixelWidth, pixelHeight);
	for (i = 0; i < rois.length; i++) {			
		roiManager("select", rois[i]);								 
		Overlay.addSelection(Roi.getStrokeColor, 2);
		Overlay.setPosition(0,0,0);
		Overlay.show();	
		if (addlabel) {
			setColor("white");
			Overlay.drawString(""+rois[i]+1, getValue("X")/pixelWidth-5, getValue("Y")/pixelHeight);			
			Overlay.add();
		}
		Overlay.show();
	}
}

function addDistOverlay(id,dist) {
	selectImage(dist);
	N = nSlices;	
	for (n=1;n<=N;n++) {
		selectImage(dist);
		setSlice(n);
		setThreshold(OUTOFBOUND+1, -OUTOFBOUND);
		run("Create Selection");
		selectImage(id);				
		run("Restore Selection");
		Overlay.addSelection("White", 1);
		Overlay.add();
		Overlay.show();
	}
}

function generateTestImage() {
	w = 1200; // width of the image
	h = 200; // height of the image
	d = 10; // size of the spots
	N = 5; // number of of roi
	n = 10; // number of spots / roi	
	D = 4*d; // cell size
	newImage("HyperStack", "32-bit composite-mode", w, h, 3, 1, 1);
	Stack.setChannel(1); run("Blue");
	Stack.setChannel(2); run("Red");
	Stack.setChannel(3); run("Green");
	//rho = newArray(0,0.125,0.25,0.5,0.75,1); // amount of colocalization for each ROI
	rho = Array.getSequence(N);
	setColor(255);
	for (roi = 0; roi < rho.length;roi++) {
		rho[roi] = rho[roi] / (rho.length-1);
		makeRectangle(roi*(w/rho.length)+d, d, (w/rho.length)-2*d, h-2*d);
		Roi.setPosition(1,1,1);		
		rx = getValue("BX");
		ry = getValue("BY");
		rw = getValue("Width");
		rh = getValue("Height");
		for (i = 0; i < n; i++) {
			x = rx + (rw-d) * random;
			y = ry + (rh-d) * random;
			t = random;
			if (t < 0.5+rho[roi]/2) {
				Stack.setChannel(2);
				makeOval(x, y, d, d);
				fill();
			}
			if (t > 0.5-rho[roi]/2) {
				Stack.setChannel(3);
				makeOval(x, y, d, d);
				fill();
			}
		}
		Stack.setChannel(1);
		makeOval(rx+rw/2-D/2, ry+rh/2-D/2, D, D);
		fill();
	}
	run("Select None");
	run("Gaussian Blur...", "sigma=2");
	run("Add Specified Noise...", "standard=5");
	rename("Test Image");
	roiManager("show all");
	run("Stack to Hyperstack...", "order=xyczt(default) channels=3 slices=1 frames=1 display=Color");
	for (i=1;i<=3;i++) {	Stack.setChannel(i); resetMinAndMax();}	
	run("Stack to Hyperstack...", "order=xyczt(default) channels=3 slices=1 frames=1 display=Composite");	
}

///////////////////////////////////////////////////////////////////////////////////
//                       Measurements
///////////////////////////////////////////////////////////////////////////////////
function createTable(name) {
	// create a table with title 'name' if it does not exist
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
		Roi.setStrokeColor("red");		
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

//////////////////////////////////////////////////////////////////////////////////
//                     Pre processing
//////////////////////////////////////////////////////////////////////////////////
function preprocess(id0, preprocess_type, scale) {
	selectImage(id0);	
	if (matches(preprocess_type, "LoG")) {		
		run("Select None");
		run("Duplicate...","title=tmp duplicate");
		run("32-bit");
		id1 = getImageID();		
		run("Square Root", "stack");
		run("Gaussian Blur...", "sigma="+scale+" stack"); 
		run("Convolve...", "text1=[-1 -1 -1\n-1 8 -1\n-1 -1 -1\n] normalize stack");
	} else if(matches(preprocess_type, "Rolling Ball")) {
		run("Select None");
		run("Duplicate...","title=tmp duplicate");
		run("32-bit");
		id1 = getImageID();
		run("Gaussian Blur...", "sigma=0.25 scaled stack");
		run("Subtract Background...", "rolling="+5*scale+" stack");
		run("Median...", "radius="+scale+" stack");
	} else {
		run("Select None");
		run("Duplicate...", "duplicate");
		id1 = getImageID();
	}
	return id1;
}

//////////////////////////////////////////////////////////////////////////////////
//                     Colocalization function
//////////////////////////////////////////////////////////////////////////////////

function getROIList() {
	/* Get the list of index of the ROI starting from 0.*/
	roi = newArray(roiManager("count"));
	roi = Array.getSequence(roiManager("count"));	
	return roi;
}

function getChannelThresholds(img, channels, logpfa) {
	/* Return an array of threshold for each channel
	 * img     : id of the window to process
	 * channel : list of channels (array)
	 * logpfa  : log probability of false alarm
	 */
	selectImage(img);		
	run("Select None"); 
	thresholds = newArray(channels.length);
	for (i = 0; i < channels.length; i++) {
		c = channels[i];		
		Stack.setChannel(c);				
		thresholds[i] = defineThreshold(logpfa);		
	}
	return thresholds;
}

function getAllchildrenROI(img, parents, channels, thresholds, areafilter, circfilter, colors) {
	/* Create ROI given the list of thresholds for each channel and return the 
	 * list of indices of ROI belonging to all ROI
	 */
	selectImage(img);
	all_children = newArray();	
	for (i = 0; i < channels.length; i++) {		
		children = addChildROI(parents, channels[i], thresholds[i], areafilter, circfilter, colors[channels[i]-1]);
		all_children = Array.concat(all_children, children);
	}
	roiManager("Show All without labels");
	return all_children;
}

function getROIChannels(rois) {
	/* returns the list of channel of the ROIs from their positions
	 * Input
	 * rois : array of ROI indices
	 * Output: array of channel for each roi
	 */
	roi_channels = newArray(rois.length);
	for (i = 0; i < rois.length; i++) {
		roiManager("select", rois[i]);
		//Roi.getPosition(channel, slice, frame);
		Stack.getPosition(channel, slice, frame);		
		roi_channels[i] = channel;
	}	
	return roi_channels;
}

function addChildROI(parents, channel, threshold, area, circ, color) {
	/* Add regions above threshold to the ROI manager and set the channel
	 * Input:
	 * parents : array of parent roi
	 * area : area filter
	 * circ : circularity filter
	 * channel : channel 
	 */
	run("Select None");
	Stack.getDisplayMode(display_mode);
	Stack.setDisplayMode("grayscale"); // to enable threshold on channel	
	// Select all parents to restrict anayluse particle to the part of the image in parents ROI	
	if (parents.length > 1) {
		roiManager("select", parents);
		roiManager("combine");
	} else {
		roiManager("select", parents[0]);
	}
	
	Stack.setChannel(channel);
	getStatistics(area, mean, min, max, std, histogram);
	setThreshold(threshold, max);
	// add new roi to the ROI manager
	n0 = roiManager("count");
	//run("Analyze Particles...", "size="+area+" pixel circularity="+circ+" add");
	run("Analyze Particles...", "add");
	n1 = roiManager("count");
	roiManager("show all");
	children = newArray(n1 - n0);
	for (i = 0; i < children.length; i++) {
		roiManager("select", n0+i);
		Roi.setPosition(channel, 1, 1);
		Roi.setStrokeColor(color);
		roiManager("update");
		children[i] = n0 + i;
	}
	Stack.setDisplayMode(display_mode);
	return children;
}

function hasROIAny(r1, r2) {
	// Test if ROI r2 belong to ROI r2
	ret = false;
	roiManager("select", r2);
	Roi.getCoordinates(x, y);
	roiManager("select", r1);	
	for (i = 0; i < x.length; i++) {
		if (Roi.contains(x[i], y[i])) {
			ret = true;
			break;
		}
	}
	return ret;
}

function hasROIAll(r1, r2) {
	// Test if ROI r2 belong to ROI r2
	ret = true;
	roiManager("select", r2);
	Roi.getCoordinates(x, y);
	roiManager("select", r1);	
	for (i = 0; i < x.length; i++) {
		if (!Roi.contains(x[i], y[i])) {
			ret = false;
			break;
		}
	}
	return ret;
}

function defineThreshold(logpfa) {
	 // Define a threshold based on a probability of false alarm (<0.5)
	 st = laplaceparam();
	 pfa = pow(10,logpfa);
	 if (pfa > 0.5) {
	 	return st[0] + st[1]*log(2*pfa);
	 } else {
	 	return st[0] - st[1]*log(2*pfa);
	 }
}

function laplaceparam() {
	// compute laplace distribution parameters 
	getHistogram(values, counts, 256);
	n = getWidth() * getHeight();
	c = 0;
	for (i = 0; i < counts.length && c <= n/2; i++) {
		c = c + counts[i];
	}
	m = values[i-1];
	b = 0;
	for (i = 0; i < counts.length; i++){
		b = b + abs(values[i] - m) * counts[i];
	}
	return newArray(m,b/n);
}

function getchildrenUnionInParents(channels, parents, children, child_channels,union_mask) {	
	/* Union of children in given parent and in either channels
	 * Output:
	 * union : the list of indices in ROI manager of selected roi
	 * union_mask : array where the indices of selected children is 1
	 */
	union = newArray(children.length);
	l = 0;
	for (i = 0; i < children.length; i++) {
		// if they belong one of the two channels and roi k
		test = false;
		for (j = 0; j < channels.length; j++) {
			test = test || (child_channels[i] == channels[j]);
		}
		if (test) {
			if (hasROIAny(parents, children[i])) {	
				union[l] = children[i];
				union_mask[i] = true;
				l++;
			}
		} else {
			union_mask[i] = false;
		}
	}
	union = Array.trim(union,l);
	return union;
}
function getIntensities() {
	// Get the intensities inside the current ROI as an array
	Roi.getContainedPoints(xpoints, ypoints);
	values = newArray(xpoints.length);
	for (i = 0; i < xpoints.length; i++) {
		values[i] = getPixel(xpoints[i], ypoints[i]);
	}
	return values;		
}


function spearman(x1,x2) {
	// https://en.wikipedia.org/wiki/Spearman%27s_rank_correlation_coefficient
	// The ranking function is not consistent with other statistical software	
	r1 = Array.rankPositions(x1);
	r2 = Array.rankPositions(x2);
	pcc = pearson(r1,r2);	
	return -pcc[0];
}

function pearson(x1, x2) {	
	// Peason correlation coefficient	
	// return an array  0:pcc,1:mu1,2:sigma1,3:mu2,4:sigma2
	if (x2.length != x1.length ) {
		print("non matching array length in prearson x1:"+x1.length+", x2:"+x2.length);
	} 
	Array.getStatistics(x1, min1, max1, mu1, sigma1);
	Array.getStatistics(x2, min1, max1, mu2, sigma2);
	n = x1.length;
	sxy = 0;
	sxx = 0;
	syy = 0;
	for (i = 0; i < x1.length; i++) {
		sxy += (x1[i]-mu1)*(x2[i]-mu2);
		sxx += (x1[i]-mu1)*(x1[i]-mu1);
		syy += (x2[i]-mu2)*(x2[i]-mu2);
	}
	pcc = sxy / (sqrt(sxx)*sqrt(syy)); 
	return newArray(pcc,mu1,sigma1,mu2,sigma2);
}

function manders(x1,x2,y1,y2,t1,t2) {
	// Manders overlap coeffcicients
	// x1: value of preprocess image in channel 1
	// x2: value of preprocess image in channel 1	
	// t1 and t2 are threshold values
	// return an array: 0:m1, 1:m2, 2:moc, 3:n1, 4:n2, 5:n12 

	if (x2.length != x1.length || y1.length != x1.length || y2.length != x1.length) {
		print("non matching array length in manders x1:"+x1.length+", x2:"+x2.length+"y1:"+y1.length+", y2:"+y2.length);		
	} 
	
	n1 = 0;
	n2 = 0;
	n12 = 0;
	s1 = 0;//sum R
	s2 = 0;//sum G
	s12 = 0;// sum Rcoloc
	s21 = 0;// sum Gcoloc
	moc12 = 0;
	moc1 = 0;
	moc2 = 0;
	for (i = 0; i < x1.length; i++) {
		if (x1[i] > t1) {n1++; s1 += y1[i]; } 
		if (x2[i] > t2) {n2++; s2 += y2[i]; } 
		if (x1[i] > t1 && x2[i] > t2) { n12++;  s12 += y1[i]; s21 += y2[i];}
		moc12 += y1[i]*y2[i];
		moc1  += y1[i]*y1[i];
		moc2  += y2[i]*y2[i];
	}
	m1 = s12 / s1;
	m2 = s21 / s2;
	moc = moc12 / sqrt(moc1 * moc2);
	return newArray(m1,m2,moc,n1,n2,n12,s1/n1,s2/n2,s12/n12,s21/n12);
}

function objectColoc(ch1, ch2, rois, roi_channels, union_mask) {	
	/* Object based colocalization
	 *  
	 */
	if (rois.length != roi_channels.length) {
		print("In objectColoc rois ("+rois.length+") and roi_channels ("+roi_channels.length+") non matching length.")
	}
	if (rois.length != union_mask.length) {
		print("In objectColoc rois ("+rois.length+") and union_mask ("+union_mask.length+") non matching length.")
	}	
	// Compute a subset for rois and roi_channels
	n = 0;
	for (i = 0; i < union_mask.length; i++) if (union_mask[i]) {
		n++;
	}
	subset = newArray(n);
	channels = newArray(n);
	j = 0;
	for (i = 0; i < union_mask.length; i++) if (union_mask[i]) {
		subset[j] = rois[i];
		channels[j] = roi_channels[i];
		j++;
	}	
	N1 = 0;
	N2 = 0;
	N12 = 0;	
	for (i = 0; i < subset.length; i++) {		
		ci = channels[i];
		// count the number of roi in each channel
		if (ci == ch1) {
			N1++;
		} else if (ci == ch2) {
			N2++;
		}
		// count the number of overlap
		for (j = i+1; j < subset.length; j++) {
			cj = channels[j];			
			if (ci!=cj && hasROIAny(subset[i], subset[j])) {
				N12++;
			}
		}
	}
	return newArray(N1,N2,N12);
}

function measureColocalization(tblname, img1, img2, channels, thresholds, parents, children, children_channel, on_original, channel_names, keys, vals) {
	/* Measure the colocalization and report by parents
	 * 
	 */
	tblname = "Colocalization.csv";
	if (!isOpen(tblname)) {
		Table.create(tblname);
	}
	
	if (children.length != children_channel.length) {
		print("Children ("+children.length+") and children_channel ("+children_channel.length+") have different size in measureColocalization");
		
	}	
	
	for (parent_id = 0; parent_id < parents.length; parent_id++) {
		
		roiManager("select", parents[parent_id]);
		roi_area = getValue("Area");
		
		for (ch1_id = 0; ch1_id < channels.length; ch1_id++) {
			
			ch1 = channels[ch1_id];
			
			for (ch2_id = ch1_id+1; ch2_id < channels.length; ch2_id++) {								
			
				ch2 = channels[ch2_id];

				union_mask = newArray(children.length);
				union = getchildrenUnionInParents(newArray(ch1,ch2), parents[parent_id], children, children_channel, union_mask);
				
				if (union.length > 0) {
				
					// measure intensity on processed image
					selectImage(img2);
					roiManager("select", union);
					roiManager("combine");
					
					Stack.setChannel(ch1);
					x1 = getIntensities();
					Stack.setChannel(ch2);
					x2 = getIntensities();
					
					// measure intensity on oringal image + min and max to detect saturation
					selectImage(img1);
					roiManager("select", union);
					roiManager("combine");
								
					Stack.setChannel(ch1);
					m1 = getValue("Max");
					y1 = getIntensities();
					Stack.setChannel(ch2);
					m2 = getValue("Max");
					y2 = getIntensities();
																		
					// Compute the colocalization measures
					if (on_original) {
						pcc = pearson(x1, x2);
						rs  = spearman(x1, x2);
						M   = manders(x1, x2,y1, y2, thresholds[ch1_id], thresholds[ch2_id]);
					} else {
						pcc = pearson(y1, y2);
						rs  = spearman(y1, y2);
						M   = manders(x1, x2, x1, x2, thresholds[ch1_id], thresholds[ch2_id]);
					}

					O = objectColoc(ch1, ch2, children, children_channel, union_mask);
								
					// saving results to table					
					row = createTable(tblname);		
					for (k = 0; k < keys.length; k++) {
						Table.set(keys[i], row, vals[i]);
					}							
					Table.set("ROI", row, parents[parent_id]+1);
					Table.set("Ch1", row, channel_names[ch1-1]);
					Table.set("Ch2", row, channel_names[ch2-1]);
					Table.set("Pearson", row, pcc[0]);
					//Table.set("Spearman", row, rs);
					Table.set("M1", row, M[0]);
					Table.set("M2", row, M[1]);
					Table.set("MOC", row, M[2]);
					Table.set("D", row, (M[5]-M[3]*M[4]/roi_area)/roi_area);

					Table.set("N1", row,  O[0]);
					Table.set("N2", row,  O[1]);
					Table.set("N12", row, O[2]);					
					Table.set("N12/N1", row, O[2] / O[0]);
					Table.set("N12/N2", row, O[2] / O[1]);
									
					Table.set("Area Object 1 [px]", row, M[3]);
					Table.set("Area Object 2 [px]", row, M[4]);
					Table.set("Area Intersection [px]", row, M[5]);
					Table.set("Intersection vs Area Object 1 [px]", row, M[3]/M[5]);
					Table.set("Intersection vs Area Object 2 [px]", row, M[4]/M[5]);
					
					Table.set("Mean Ch1 in Object 1", row, M[6]);
					Table.set("Mean Ch2 in Object 2", row, M[7]);
					Table.set("Mean Ch1 in Intersection", row, M[8]);
					Table.set("Mean Ch2 in Intersection", row, M[9]);
					
					Table.set("Mean Ch1", row, pcc[1]);
					Table.set("Stddev Ch1 ", row, pcc[2]);
					Table.set("Max Ch1", row, m1);
					
					Table.set("Mean Ch2", row, pcc[3]);
					Table.set("Stddev Ch2", row, pcc[4]);
					Table.set("Max Ch2", row, m2);					
					
				}	else {
					print("empty set");
				}	
			}
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


function getDistROI(dist) {
	/*Get the roi corresponding to each slice of distance function dist*/
	selectImage(dist);
	N = nSlices;
	roi = newArray(N);
	run("Select None");
	for (n = 1; n <= N; n++) {
		selectImage(dist);
		setSlice(n);
		setThreshold(OUTOFBOUND+1, -OUTOFBOUND);
		run("Create Selection");
		if (0) {
			roiManager("select", n-1);
			run("Restore Selection");
			roiManager("update");
			roi[n-1] = n-1;
		} else {
			Roi.setStrokeColor("white");
			roiManager("add");			
			roi[n-1] = roiManager("count")-1;		
		}
		resetThreshold;
	}
	run("Select None");
	return roi;
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

function segmentBackground(id, channel, zoom) {
	/* Segment the background and return the image ID */
	print(" - Segment background on channel " + channel);
	selectImage(id);
	if (matches(channel, "max")) {
		run("Z Project...", "projection=[Max Intensity]");
	} else {
		run("Duplicate...", "duplicate channels="+channel);		
	}			
	run("Median...", "radius="+Math.ceil(1+5/zoom));	
	setImageThreshold(0.5);
	run("Convert to Mask");	
	run("Fill Holes");
	return getImageID;
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

function defineThreshold(logpfa) {
	 // Define a threshold based on a probability of false alarm (<0.5)
	 st = laplaceparam();
	 pfa = pow(10,logpfa);
	 if (pfa > 0.5) {
	 	return st[0] + st[1]*log(2*pfa);
	 } else {
	 	return st[0] - st[1]*log(2*pfa);
	 }
}

function laplaceparam() {
	// compute laplace distribution parameters 
	getHistogram(values, counts, 256);
	n = getWidth() * getHeight();
	c = 0;
	for (i = 0; i < counts.length && c <= n/2; i++) {
		c = c + counts[i];
	}
	m = values[i-1];
	b = 0;
	for (i = 0; i < counts.length; i++){
		b = b + abs(values[i] - m) * counts[i];
	}
	return newArray(m,b/n);
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


