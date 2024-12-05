#@File(label="Image file", description="input file (use batch for several files)") image_file
#@Integer(label="Bandwidth",value=21, description="Window size used to estimate baseline level (F0)") bandwidth
#@String(label="Metric", choices={"signal","sbr","zscore"}) metric
#@Boolean(label="Prefilter", value=true, description="Prefilter the signal") prefilter
#@Float(label="Threshold", value=0.025) threshold
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
 
function createTest() {
	run("Close All");
	if (isOpen("ROI Manager")) { close("ROI Manager"); }
	n = 100;
	m = 3;
	newImage("HyperStack", "16-bit grayscale-mode", n, n, 1, 1, 100);
	for (k = 0; k < m; k++) {
		makeOval(n*(k+1)/(m+1)-5, n/2-5, 10, 10);
		roiManager("add");
	}
	values = newArray(roiManager("count"));
	for (frame = 0; frame < 100; frame++) {
		for (i = 0; i < roiManager("count"); i++) {
			if (random < 0.1 && values[i] < 1000/4) {
				values[i] = 1000;
			}
			values[i] = values[i] * 0.5;
			roiManager("select", i);
			Stack.setFrame(frame);
			setColor(values[i]);
			fill();
			run("Select None");
		}
	}
	run("Macro...", "code=v=100+v*exp(-z/500) stack");
	
	run("Add Specified Noise...", "stack standard=5");
	run("Properties...", "channels=1 slices=1 frames=100 pixel_width=0.1 pixel_height=0.1 voxel_depth=0.1 frame=[1 sec]");
}

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


function filterAdaptive(f,sigma,niter) {
	/* 
	 *  Smooth a 1d signal
	 *  
	 *  Args:
	 *   f (array) : original signal
	 *   sigma (number) : noise standard deviation
	 *   niter (number) : number of iterations
	 *   
	 * Returns (array) filtered signal
	 */
	n = f.length;
	u0 = Array.copy(f);
	u1 = newArray(n);
	v0 = newArray(n);
	v1 = newArray(n);
	stop = newArray(n);

	for (i = 0; i < n; i++) {
		v0[i] = sigma * sigma;
	}

	for (iter = 0; iter < niter; iter++) {

		h = pow(2,iter);

		for (i = 0; i < n; i++) if (stop[i]==0) {

			j0 = maxOf(0, i - h);
			j1 = minOf(n - 1, i + h);

			swu = 0;
			swv = 0;
			sw = 0;

			for (j = j0; j <= j1; j++) {

				d = (u0[j] - u0[i]) * (u0[j] - u0[i]) /(9 * v0[i]);
				w = exp(-0.5 * d);
				sw += w;
				swu += w * f[j];
				swv += w *w *sigma*sigma;
			}

			if (sw > 0) {
				u1[i] = swu / sw;
				v1[i] = swv / sw;
			}
			
			if (abs(u1[i] - u0[i]) >   sqrt(v0[i])) {
				u1[i] = u0[i];
				v1[i] = v0[i];
				stop[i] = iter;
			}
		}

		u0 = Array.copy(u1);
		v0 = Array.copy(v1);
	}
	return u1;
}

function median(x) {
	/* Median
	 *  x (array): list of values
	 * Return the median of the array x
	 */
	y = Array.copy(x);
	Array.sort(y);
	n = y.length;
	if (n==0) {
		return 0;
	}
	if (n%2==0) {
		return 0.5*(y[n/2-1] + y[n/2]);
	} else {
		return y[(n-1)/2];
	}
}

function estimateNoiseStd(f) {
	/* Estimate noise standard deviation
	 *  f (array) : 1d signal
	 * Returns the estimated noise standard deviation
	 */
	n = f.length;
	delta = newArray(n-2);	
	for (i = 1; i < n-1; i++) {
		delta[i-1] = minOf(f[i+1] - f[i], f[i] - f[i-1]);	
	}
	m = median(delta);
	for (i = 0; i < delta.length; i++) {
		delta[i] = abs(delta[i] - m);
	}
	m = 1.48 * median(delta) / sqrt(2);
	return m;
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
	if (stdDev > 1e-6) {
		for (i = 0; i < data.length; i++) {
			dest[i] = (data[i] - mean) / stdDev;
		}
	} else {
		for (i = 0; i < data.length; i++) {
			dest[i] = (data[i] - mean);
		}
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

function findPeakEvents(score, baseline, threshold) {
	/*
	 * Find events as peaks
	 *
	 * Args:
	 *  score (array)
	 *  tolerance (float)
	 *  
	 * Returns:
	 * 	Array with value equal to maximum of peak in original signal and 0 otherwise
	 */
	events = newArray(score.length);
	for (j = 0; j < score.length-1; j++) {
		if (score[j] > threshold && score[j] >= baseline[j]) {
			events[j] = score[j];
		}
	}
	return events;
}

function findOnEvents(score, baseline, threshold) {
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
	start = -1; // start of event
	peak  = 0; // peak flag
	duration = 0; // duration of event
	amplitude = 0; // amplitude
	
	for (j = 0; j < score.length-1; j++) {
		
		// start of an event
		if (score[j] > threshold && start == -1) {
			start = j;
		}
		
		// during the event
		if (score[j] > threshold && start > -1) {
			duration++;
			if (score[j] >= baseline[j]) {
				peak++;
				if (score[j] > amplitude && peak == 1) {
					amplitude = score[j];
				}
			}
		}
		
		// end of events: under baseline or peak > 1
		if ((score[j] < threshold || peak > 1 || j >= score.length-2) && start > -1 ) {
			if (peak > 1 && duration > 2) {
				duration = duration - 3;
				j = j - 1;
			}
			// set the event array to the max amplitude
			for (k = 0; k <= duration; k++) {
				events[start+k] = amplitude ;
			}
			// reset the events accumulators
			start = -1;
			duration = 0;
			amplitude = 0;
			peak = 0;
		}
	}
	return events;
}

function computeAllSignals(table_name, bandwidth, metric, prefilter, tolerance) {
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
	 		f = Array.copy(zscore);
	 	} else if (metric == "sbr") {
	 		f = Array.copy(sbr);
	 	} else {
	 		f = Array.copy(signal);
	 	}
	 	if (prefilter) {
	 		sigma = estimateNoiseStd(f);
	 		f = filterAdaptive(f, sigma, 4);
	 	}
	 	baseline = filterRank(f, 5, 1);
	 	peaks = findPeakEvents(f, baseline, tolerance);
	 	onoff = findOnEvents(f, baseline, tolerance);
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
			if ((onoff[j] == 0 || j == onoff.length - 1) && start > -1) {
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
		n = 0; // number of frames
		s = 0; // sum of amplitude
		d = 0; // sum of durations
		level = 0;
		duration = 0;
		for (j = 0; j < onoff.length-1; j++) {
			// start segment
			if (onoff[j] > 0 && level == 0) {
				n += 1;
				level = 1;
			}
			// segment
			if (onoff[j] > 0 && level == 1) {
				duration++;
			}
			// end of segment
			if (onoff[j] != onoff[j+1] && level == 1) {
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
	w = 680;
	h = 400;
	newImage("Figure", "RGB white", w, h, nroi);
	dst = getImageID();
	// compute min max 
	yMin = 1e6;
	yMax = 0;
	for (roi = 0; roi < nroi; roi++) {
		selectWindow(tbl1);
		if (metric == "signal") {
			sig = Table.getColumn(columns[K*roi+1]);
		} else if (metric == "sbr") {
			sig = Table.getColumn(columns[K*roi+2]);
		} else {
			sig = Table.getColumn(columns[K*roi+3]);
		}
		Array.getStatistics(sig, min, max, mean, stdDev);
		yMin = minOf(yMin, min);
		yMax = maxOf(yMax, max);
	}
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
		Plot.setLimits(0, sig.length, yMin, yMax);
		Plot.update();
		selectWindow("plot");
		makeRectangle(20, 0, w, h);
		run("Copy");
		selectImage(dst);
		Stack.setSlice(roi+1);
		//makeRectangle(0, roi*h, w, h);
		run("Paste");
		selectWindow("plot");close();
	}
	selectImage(dst);
	run("Make Montage...", "columns="+round(sqrt(nroi))+" rows="+Math.ceil(nroi/round(sqrt(nroi)))+" scale=1");
	mnt = getImageID();
	selectImage(dst);close();
	return mnt;
}

function segmentActiveROI() {
	run("Z Project...", "projection=[Max Intensity]");
	run("Median...", "radius=5");
	setThreshold(getValue("Mean")+3*getValue("StdDev"), getValue("Max"));
	run("Analyze Particles...", "size=100-Infinity clear add");
	close();
}


// close all images
run("Close All");
var batchid;
batchid = batchid + 1;
print("~~~~~~~~~~["+batchid+"]~~~~~~~~~~~~~~~");

if (matches(image_file, ".*test")>0) {
	createTest();
	interval = 1;
	frames = 100;
	unit = "s";
	rename("test");
	is_test = true;
} else {
	is_test = false;
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
				//print("Segmenting active ROIs");
				//segmentActiveROI();
				print("Use the whole image");
				run("Select All");
				roiManager("Add");
				if (roiManager("count")==0) {
					print("ERROR: No ROI");
					exit();
				}
			}
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
computeAllSignals(table_signal, bandwidth, metric, prefilter, threshold);

fig = makeFigure(table_signal, metric);

// Store events properties in a table
reportPeakEvents(table_signal, table_peak_events, interval, frames, unit);
reportOnEvents(table_signal, table_onoff_events, interval, frames, unit);

// report the events statistics in a table
computePeakEventsStats(table_signal, table_peak_stats, interval, frames, unit);
computeOnEventsStats(table_signal, table_onoff_stats, interval, frames, unit);

// Save the 3 tables
if (save_tables && !is_test) {
	selectWindow(table_signal);
	print("Saving signals:\n" +  folder + File.separator +table_signal);
	Table.save(folder + File.separator + table_signal);
	run("Close");
	
	selectWindow(table_peak_events);
	print("Saving events:\n" +  folder + File.separator + table_peak_events);
	Table.save(folder + File.separator + table_peak_events);
	run("Close");
	
	selectWindow(table_peak_stats);
	print("Saving statistics:\n" +  folder + File.separator + table_peak_stats);
	Table.save(folder + File.separator + table_peak_stats);
	run("Close");
	
	selectWindow(table_onoff_events);
	print("Saving events:\n" +  folder + File.separator + table_onoff_events);
	Table.save(folder + File.separator + table_onoff_events);
	run("Close");
	
	selectWindow(table_onoff_stats);
	print("Saving statistics:\n" +  folder + File.separator + table_onoff_stats);
	Table.save(folder + File.separator + table_onoff_stats);
	run("Close");
	
	figname = folder + File.separator + File.getNameWithoutExtension(image_file) + "-figure.jpg";
	print("Saving figure\n"+ figname);
	selectImage(fig);
	saveAs("JPG", figname);
	close();
}
run("Select None");
setBatchMode("exit and display");
