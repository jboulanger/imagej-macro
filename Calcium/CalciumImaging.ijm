#@File(label="image file") image_file
#@Integer(label="bandwidth",value=21) bandwidth
#@String(label="metric", choices={"signal","sbr","zscore"}) metric
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
 * Jerome for Feline (Sep-2024), Melissa (Nov-2024)
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

function findPeakEvents(score, threshold) {
	/*
	 * Find events as peaks
	 *
	 * Args:
	 *  score (array)
	 *  tolerance (float)
	 * Returns:
	 * 	Array with value equal to amplitude of peak and 0 otherwise
	 */
	events = newArray(score.length);
	smoothed = filterRank(score, 5, 0.5);
	baseline = filterRank(score, 7, 1); // maximum filter
	for (j = 1; j < score.length-1; j++) {
		if (smoothed[j] > threshold && score[j] >= baseline[j-1] && score[j] >= baseline[j+1]) {
			events[j] = smoothed[j];
		}
	}
	return events;
}

function findOnEvents(score, threshold) {
	/*
	 * Identify events as ON/OFF events
	 *
	 * Peak needs to be separated
	 *
	 * Args:
	 * 	score (array)
	 * 	tolerance (float)
	 * Returns:
	 * 	Array with mean amplitude and 0 otherwise
	 */
	events = newArray(score.length);
    //baseline = filterRank(score, 7, 1);
	start = -1; // start of event
	duration = 0; // duration of event
	amplitude = 0; // amplitude
	for (j = 0; j < score.length-1; j++) {
		// start of an event
		if (score[j] > threshold && start == -1) {
			start = j;
		}
		if (score[j] > threshold && start > -1) {
			amplitude += score[j];
			duration++;
		}
		// end of events
		if (score[j] < threshold && start > -1) {
			// set the event array to the average amplitude
			for (k = 0; k <= duration; k++) {
				events[start+k] = amplitude / duration;
			}
			// reset the events accumulators
			start = -1;
			duration = 0;
			amplitude = 0;
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
	 		peaks = findPeakEvents(zscore, tolerance);
	 		onoff = findOnEvents(zscore, tolerance);
	 	} else if (metric == "sbr") {
	 		peaks = findPeakEvents(sbr, tolerance);
	 		onoff = findOnEvents(sbr, tolerance);
	 	} else {
	 		peaks = findPeakEvents(signal, tolerance);
	 		onoff = findOnEvents(signal, tolerance);
	 	}
	 	Table.setColumn(name+"-signal", signal);
	 	Table.setColumn(name+"-sbr", sbr);
	 	Table.setColumn(name+"-zscore", zscore);
	 	Table.setColumn(name+"-peaks", peaks);
	 	Table.setColumn(name+"-on", onoff);
	 	Table.update;
	 }
	 return 0;
}

function reportPeakEvents(tbl1, tbl2, interval, frames, unit) {
	/*
	 * Store information for each peak event
	 */
	selectWindow(tbl1);
	columns = split(Table.headings,"\t");
	Table.create(tbl2);
	K = 5;
	for (i = K-1; i < columns.length; i+=K) {
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

function reportOnEvents(tbl1, tbl2, interval, frames, unit) {
	/*
	 * Store information for each ON/OFF events
	 * 
	 * Args: 
	 *	tbl1 (str) : signal table
	 *	tbl2 (str) : peak eavent table
	 *	interval (float): time interval
	 *	frames (int): number of frames
	 *	unit (str): unit of time
	 */
	K = 5;
	selectWindow(tbl1);
	columns = split(Table.headings,"\t");
	frame = Table.getColumn("Frame");
	nroi = (columns.length - 1) / K;
	Table.create(tbl2);
	for (roi = 0; roi < nroi; roi++) {
		selectWindow(tbl1);
		onoff = Table.getColumn(columns[K*roi+5]);
		peaks = Table.getColumn(columns[K*roi+4]);
		name = replace(columns[K*roi+4], "-peaks", "");
		event_id = 0;
		start = -1;
		duration = 0;
		npeaks = 0;
		for (j = 0; j < onoff.length; j++) {
			// start
			if (onoff[j] > 0 && start == -1) {
				start = j;
				duration = 0;
			} 
			// segment
			if (onoff[j] > 0 && start > -1) {
				duration++;
				npeaks += peaks[j];
			}
			// stop
			if (onoff[j] == 0 && start > -1) {
				selectWindow(tbl2);
				row = Table.size;
				Table.set("roi", row, name);
				Table.set("event", row, event_id);
				Table.set("frame start", row, start);
				Table.set("frame end", row, j-1);
				Table.set("time (s)", row, start * interval);
				Table.set("duration (s)", row, duration * interval);
				Table.set("mean amplitude", row, onoff[j-1]);
				event_id++;
				start = -1;
				duration = 0;
				npeaks = 0;
			}
		}
		Table.update;
	}
}


function computePeakEventsStats(tbl1, tbl2, interval, frames, unit) {
	/*
	 * Compute the number of peaks, frequency and mean amplitude
	 * 
	 * Args: 
	 *	tbl1 (str) : signal table
	 *	tbl2 (str) : peak eavent table
	 *	interval (float): time interval
	 *	frames (int): number of frames
	 *	unit (str): unit of time
	 */
	K = 5;
	selectWindow(tbl1);
	columns = split(Table.headings,"\t");
	nroi = (columns.length - 1) / K;
	Table.create(tbl2);
	for (roi = 0; roi < nroi; roi++) {
		selectWindow(tbl1);
		peaks = Table.getColumn(columns[K*roi+4]);
		// count the number of peaks
		n = 0;
		s = 0;
		for (j = 0; j < peaks.length; j++) if (peaks[j] > 0) {
			n += 1;
			s += peaks[j];
		}
		roiManager("select", roi);
		selectWindow(tbl2);
		name = replace(columns[K*roi+4], "-peaks", "");
		row = Table.size;
		Table.set("roi", row, name);
		Table.set("x",row, getValue("X"));
		Table.set("y",row, getValue("Y"));
		Table.set("number", row, n);
		Table.set("frequency [mHz]", row, 1000*n / (interval * frames));
		Table.set("mean amplitude", row, s/n);
		Table.update;
	}
}


function computeOnEventsStats(tbl1, tbl2, interval, frames, unit) {
	/*
	 * Compute the number of peaks, frequency and mean amplitude per roi
	 * 
	 * Args: 
	 *	tbl1 (str) : signal table
	 *	tbl2 (str) : on/off event table
	 *	interval (float): time interval
	 *	frames (int): number of frames
	 *	unit (str): unit of time
	 */
	K = 5;
	selectWindow(tbl1);
	columns = split(Table.headings,"\t");
	nroi = (columns.length - 1) / K;
	Table.create(tbl2);
	for (roi = 0; roi < nroi; roi++) {
		selectWindow(tbl1);
		onoff = Table.getColumn(columns[K*roi+5]);
		name = replace(columns[K*roi+5], "-on", "");
		// count the leading front of peaks
		n = 0;
		s = 0;
		d = 0;
		level = 0;
		duration = 0;
		for (j = 0; j < onoff.length-1; j++) {
			// start segment
			if (onoff[j] > 0 && level == 0) {
				n += 1;
				s += onoff[j];
				level = 1;
			}
			// segment
			if (onoff[j] > 0 && level == 1) {
				duration++;
			}
			// end of segment
			if (onoff[j] == 0 && level == 1) {
				d += duration;
				s += onoff[j-1];
				level = 0;
				duration = 0;
			}
		}
		selectWindow(tbl2);
		roiManager("select", roi);
		row = Table.size;
		Table.set("roi", row, name);
		Table.set("x",row, getValue("X"));
		Table.set("y",row, getValue("Y"));
		Table.set("number", row, n);
		Table.set("frequency [mHz]", row, 1000*n / (interval * frames));
		Table.set("mean amplitude", row, s/n);
		Table.set("mean duration", row, d/n * interval);
		Table.update;
	}
}

function makeFigure(tbl1, metric) {
	/* 
	 *  Make a figure
	 *  
	 *  Args:
	 *   tbl1 (str): name of the table
	 *   metric (str): name of the metric for the graphs
	 *  
	 */
	K = 5;
	selectWindow(tbl1);
	columns = split(Table.headings,"\t");
	frame = Table.getColumn("Frame");
	nroi = (columns.length - 1) / K;
	print("nroi:", nroi);
	w = 680;
	h = 400;
	newImage("Figure", "RGB white", w, h*nroi, 1);
	dst = getImageID();
	for (roi = 0; roi < nroi; roi++) {
		selectWindow(tbl1);
		if (metric == "signal") {
			sig = Table.getColumn(columns[K*roi+1]);
		} else if (metric == "sbr") {
			sig = Table.getColumn(columns[K*roi+2]);
		} else {
			sig = Table.getColumn(columns[K*roi+3]);
		}
		peaks = Table.getColumn(columns[K*roi+4]);
		onoff = Table.getColumn(columns[K*roi+5]);
		Plot.create("plot", "Time", metric);
		Plot.setColor("red");
		Plot.add("line", frame, peaks);
		Plot.setColor("blue");
		Plot.add("line", frame, sig);
		Plot.setColor("green");
		Plot.add("line", frame, onoff);
		Plot.setLimitsToFit();
		Plot.update();
		selectWindow("plot");
		makeRectangle(20, 0, w, h);
		run("Copy");
		selectImage(dst);
		makeRectangle(0, roi*h, w, h);
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

// bioformat does not decode the frame interval of hyperstacks
if (endsWith(image_file, ".tif") || endsWith(image_file, ".tiff")) {
	open(image_file);
} else {
	run("Bio-Formats Importer", "open=["+image_file+"] color_mode=Composite rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
}

// get the frame interval, unit and number
selectWindow(File.getName(image_file));
interval = Stack.getFrameInterval();
Stack.getUnits(X, Y, Z, unit, Value);
Stack.getDimensions(width, height, channels, slices, frames);
print("Image width:",width,", height:",height, ", depth:", slices, ", frames:", frames);

if (frames==1 && slices > 1) {
	print("WARNING: the number of frame is 1. Use Stack to Hyperstack to fix the dimensions.");
}

if (interval == 0) {
	print("WARNING: frame interval is 0, will not be able to compute frequencies.");
}

// Close the ROI manager
if (isOpen("ROI Manager")) {
	selectWindow("ROI Manager");
	run("Close");
}

// find the roi file
roi_file = folder + File.separator + File.getNameWithoutExtension(image_file) + "_ROIset.zip";
if (File.exists(roi_file)) {
	print("Opening " + File.getName(roi_file));
	roiManager("open", roi_file);
} else {
	roi_file = folder + File.separator + File.getNameWithoutExtension(image_file) + "_ROIset.roi";
	if (File.exists(roi_file)) {
		print("Opening " + File.getName(roi_file));
		roiManager("open", roi_file);
	} else {
		if (roiManager("count")==0) {
			print("ERROR: No ROI");
			exit();
		}
	}
} 

print("Number of ROIs", roiManager("count"));


// start processing
setBatchMode("hide");

// define the name of the output tables
table_signal = File.getNameWithoutExtension(image_file) + "-signal.csv";
table_peak_events = File.getNameWithoutExtension(image_file) + "-peak-events.csv";
table_peak_stats = File.getNameWithoutExtension(image_file) + "-peak-stats.csv";
table_onoff_events = File.getNameWithoutExtension(image_file) + "-onoff-events.csv";
table_onoff_stats = File.getNameWithoutExtension(image_file) + "-onoff-stats.csv";

// process all ROI
computeAllSignals(table_signal, bandwidth, metric, threshold);

fig = makeFigure(table_signal, metric);

// Store events properties in a table
reportPeakEvents(table_signal, table_peak_events, interval, frames, unit);
reportOnEvents(table_signal, table_onoff_events, interval, frames, unit);

// report the events statistics in a table
computePeakEventsStats(table_signal, table_peak_stats, interval, frames, unit);
computeOnEventsStats(table_signal, table_onoff_stats, interval, frames, unit);

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
run("Select None");
setBatchMode("exit and display");
