#@File(label="image file") image_file
#@Integer(label="bandwidth",value=21) bandwidth
#@String(label="metric", choices={"sbr","zscore"}) metric
#@Float(label="threshold", value=0.025) threshold
#@Boolean(label="Save Tables", value=true) save_tables
#@OUTPUT number_of_events

/*
 * Analysis of calcium signal 
 * 
 * Compute SBR as (F-F0)/F0
 * where F0 is the baseline estimated as a 1d median filter 
 * followed by a 1d minimum filter.
 * 
 * Events are detected as SBR or zscore local maxima above a threshold.
 * 
 * Jerome for Feline (Sep-2024)
 */

function computeSignal() {
	/*
	 * Extract time signal for current image and region
	 * 
	 * Returns
	 *  array with mean intensity over time
	 */
	Stack.getDimensions(width, height, channels, slices, frames);
	values = newArray(frames);
	for (frame = 1; frame <= frames; frame++) {
		Stack.setFrame(frame);
		values[frame-1] = getValue("Mean");
	}
	return values;
}

function correctBleaching() {
	/* Correct bleaching using an exponential model	*/
	run("Select None");
	run("32-bit");
	values = extractSignal();
	time = Array.getSequence(values.length);
	Fit.doFit("Exponential recovery", time, values);	
	b = Fit.p(1);
	c = Fit.p(2);
	Stack.getDimensions(width, height, channels, slices, frames);
	run("Macro...","v=(v-("+c+"))*(exp(("+b+")*z)) stack");		
}

function filterRank(values, bandwidth, rank) {
	/*
	 * Filter the signal with a rank filter
	 * Args:
	 *  values (array): signal to process
	 *  bandwidth (int): bandwidth of the filter
	 *  rank (float): rank
	 * Returns
	 *  array: filtered data
	*/
	n = values.length;
	dest = newArray(n);
	for (i = 0; i < n; i++) {			
		k_start = maxOf(0, i - Math.round((bandwidth-1)/2));
		k_end = minOf(n, i + Math.round((bandwidth-1)/2) + 1);		
		segment = Array.slice(values, k_start, k_end);
		Array.sort(segment);
		dest[i] = segment[maxOf(0,minOf(segment.length-1,Math.floor(rank*segment.length)))];
	}
	return dest;
}

function subtractArray(x, y) {
	/* 
	 * Subtract two arrays (x-y)
	 * 
	 * Args:
	 *  x (array) 
	 *  y (array)
	 * Returns
	 *  z (array)
	 */
	if (x.length != y.length) {
		exit("array1 and array2 have different length");
	}
	n = x.length;
	z = newArray(n);
	for (i = 0; i < n; i++) {
		z[i] = x[i] - y[i];
	}
	return z;
}

function ratioArray(x, y) {
	/*
	 * 	Ratio of two arrays
	 * 	
	 * Args:
	 *  x (array) 
	 *  y (array)
	 * Returns
	 *  z (array)
	 * 		 
	 */
	if (x.length != y.length) {
		exit("x and y have different length");
	}
	n = x.length;
	z = newArray(n);
	for (i = 0; i < n; i++) {
		z[i] = x[i] / y[i];
	}
	return z;
}

function standardizeArray(data) {
	/*
	 * Standadize the array (data-mean)/std
	 * 
	 * Args:
	 *   data (array): input data to normalize
	 * Returns:
	 *   array: standardized array
   	 */
	Array.getStatistics(data, min, max, mean, stdDev);
	dest = newArray(n);
	for (i = 0; i < data.length; i++) {
		dest[i] = (data[i] - mean) / stdDev;
	}
	return dest;
}

function computeSBR(data, bandwidth) {
	/*
	 * Compute signal to base ratio (SBR)
	 *
	 */
	 
	// apply a median filter
	med = filterRank(data, bandwidth, 0.5);
	
	// apply a minumum filter
	F0 = filterRank(med, bandwidth, 0);
	
	// compute the ratio (F-F0)/F0
	ratio = ratioArray(subtractArray(data,  F0), F0);
	return ratio;
}

function findEvents(score, threshold) {
	/*
	 * Find events
	 * 
	 * Args
	 *  score (array)
	 *  tolerance (float)
	 * Returns
	 * 	Array
	 */
	events = newArray(score.length);
	baseline = filterRank(score, 7, 1);
	for (j = 1; j < score.length-1; j++) {		
		if (score[j] > threshold && score[j] >= baseline[j-1] && score[j] >= baseline[j+1]) {
			events[j] = score[j];
		}
	}
	return events;
}

function computeAllSignals(table_name, bandwidth, metric, tolerance) {	
	/* 
	 * Extract signal from image for each ROI
	 * and store them in a Table
	 */
	 Table.create(table_name);
	 Stack.getDimensions(width, height, channels, slices, frames);
	 time = Array.getSequence(frames);
	 Table.setColumn("Frame", time);
	 n = roiManager("count");
	 for (i = 0; i < n; i++) {
	 	roiManager("select", i);
	 	name = Roi.getName();
	 	print("Extracting intensity from ROI ", name);
	 	signal = computeSignal();
	 	sbr = computeSBR(signal, bandwidth);
	 	zscore = standardizeArray(sbr);
	 	if (metric == "zscore") {
	 		peaks = findEvents(zscore, tolerance);	 	
	 	} else {
	 		peaks = findEvents(sbr, tolerance);	 	
	 	}
	 	Table.setColumn(name+"-signal", signal);
	 	Table.setColumn(name+"-sbr", sbr);
	 	Table.setColumn(name+"-zscore", zscore);	 	
	 	Table.setColumn(name+"-peaks", peaks);
	 	Table.update;
	 }
	 return 0;
}

function reportEvents(tbl1, tbl2, interval, frames, unit) {
	/*
	 * Store information for each event
	 */
	selectWindow(tbl1);
	columns = split(Table.headings,"\t");
	Table.create(tbl2);
	K = 4;	
	for (i = K; i < columns.length; i+=K) {
		selectWindow(tbl1);		
		peaks = Table.getColumn(columns[i]);				
		name = replace(columns[i], "-peaks", "");
		event_id = 0;
		for (j = 0; j < peaks.length; j++) if (peaks[j] > 0) {
			selectWindow(tbl2);
			row = Table.size;
			Table.set("roi", row, name);
			Table.set("event", row, event_id);
			Table.set("frame", row, j);
			Table.set("time", row, j * interval);
			Table.set("amplitude", row, peaks[j]);
			event_id++;
		}		
		Table.update;
	}
}
	

function computeEventsStats(tbl1, tbl2, interval, frames, unit) {
	/* 
	 * Compute the number of peaks 
	 * 
	 */
	selectWindow(tbl1);
	columns = split(Table.headings,"\t");
	Table.create(tbl2);
	K = 4;	
	for (i = K; i < columns.length; i+=K) {
		selectWindow(tbl1);		
		peaks = Table.getColumn(columns[i]);		
		n = 0;
		s = 0;
		for (j = 0; j < peaks.length; j++) if (peaks[j] > 0) {			
			n += 1;
			s += peaks[j];
		}	
		selectWindow(tbl2);
		name = replace(columns[i], "-peaks", "");
		row = Table.size;
		Table.set("roi", row, name);
		Table.set("number", row, n);
		Table.set("frequency [mHz]", row, 1000*n / (interval * frames));		
		Table.set("mean amplitude", row, s/n);
		Table.update;
	}
}

function makeFigure(tbl1) {
	/*Make a figure*/	
	K = 4;	
	selectWindow(tbl1);
	columns = split(Table.headings,"\t");	
	frame = Table.getColumn("Frame");
	nroi = (columns.length - 1) / 4;
	w = 680;
	h = 400;
	newImage("Figure", "RGB white", w, h*nroi, 1);
	dst = getImageID();
	for (i = K; i < columns.length; i+=K) {				
		selectWindow(tbl1);		
		peaks = Table.getColumn(columns[i]);
		sbr = Table.getColumn(columns[i-2]);		
		Plot.create("plot", "Time", "SBR");
		Plot.setColor("red");
		Plot.add("line", frame, peaks);
		Plot.setColor("blue");
		Plot.add("line", frame, sbr);			
		Plot.setLimitsToFit();
		Plot.update();
		selectWindow("plot");
		makeRectangle(20, 0, w, h);
		run("Copy");	
		selectImage(dst);
		makeRectangle(0, (i-K)/K*h, w, h);
		run("Paste");
		selectWindow("plot");close();
	}
	return dst;
}

// close all images
run("Close All");
var batchid;
batchid = batchid + 1;
print("~~~~~~~~~~["+batchid+"]~~~~~~~~~~~~~~~");
folder = File.getDirectory(image_file);
print("Current folder :" + folder);

// open the image
print("Opening " + File.getName(image_file));
run("Bio-Formats Importer", "open=["+image_file+"] color_mode=Composite rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");

// get the frame interval, unit and number
interval = Stack.getFrameInterval();
Stack.getUnits(X, Y, Z, unit, Value);
Stack.getDimensions(width, height, channels, slices, frames);

// Close the ROI manager
if (isOpen("ROI Manager")) {
	selectWindow("ROI Manager");
	run("Close");
}

// find the roi file automatically
roi_file = folder + File.separator + File.getNameWithoutExtension(image_file) + "_ROIset.zip";
if (!File.exists(roi_file)) {		
	roi_file = folder + File.separator + File.getNameWithoutExtension(image_file) + "_ROIset.roi";	
} 
if (!File.exists(roi_file)) {		
	print("ROI file not found " + roi_file);
	exit();
}
print("Opening " + File.getName(roi_file));
open(roi_file);
print("Number of ROIs", roiManager("count"));

// start processing
setBatchMode("hide");

// define the name of the output tables
table_signal = File.getNameWithoutExtension(image_file) + "-signal.csv";
table_events = File.getNameWithoutExtension(image_file) + "-events.csv";
table_stats = File.getNameWithoutExtension(image_file) + "-stats.csv";

// process all ROI
computeAllSignals(table_signal, bandwidth, metric, threshold);

fig = makeFigure(table_signal);

// Store events properties in a table
reportEvents(table_signal, table_events, interval, frames, unit);

// report the events statistics
computeEventsStats(table_signal, table_stats, interval, frames, unit);

// Save the 3 tables
if (save_tables) {
	selectWindow(table_signal);
	print("Saving signals:\n" +  folder + File.separator +table_signal);
	Table.save(folder + File.separator + table_signal);
	run("Close");
	selectWindow(table_events);
	print("Saving events:\n" +  folder + File.separator + table_events);
	Table.save(folder + File.separator + table_events);
	run("Close");
	selectWindow(table_stats);
	print("Saving statistics:\n" +  folder + File.separator + table_stats);
	Table.save(folder + File.separator + table_stats);
	run("Close");
	figname = folder + File.separator + File.getNameWithoutExtension(image_file) + "-figure.jpg";
	print("Saving figure\n"+ figname);	
	selectImage(fig);
	saveAs("JPG", figname);
	close();
}

setBatchMode("exit and display");
