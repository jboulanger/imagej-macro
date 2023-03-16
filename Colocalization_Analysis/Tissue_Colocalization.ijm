//@String(label="channel names",value="A,B,C,D", description="comma separated list of names") channel_names_str 
//@Integer(label="Channel with nuclei labeling for segmentation", value=1) channel_nuclei

channel_nuclei = 1;
run("Close All");
open("/home/jeromeb/Desktop/smallcrop.tif");

makePolygon(480,57,84,53,80,440,349,435,488,278,480,154);

filename = getTitle();
tbl1 = "Measure per ROI.csv";
tbl2 = "Summary.csv";
closeWindow(tbl1);
closeWindow("ROI Manager");
id = getImageID();

setBatchMode("hide");
channel_names = parseCSVString(channel_names_str);
cropActiveSelection();
rois = segmentNuclei(id, channel_nuclei);

measureInROI(tbl1, id, rois, channel_names);

masks = mapPostiveROIs(id, rois);

recordPositiveROIs(tbl1, masks, rois);

measureROIdistanceToLabels(tbl1, masks, rois);

summarizeTable(tbl1, tbl2, filename, channel_names);

selectImage(masks); close();

selectImage(id);
run("Select None");
addROIsToOverlay(id, rois);
run("Collect Garbage");	
run("Select None");
setBatchMode("exit and display");
print("Done");
 

function segmentNuclei(id, channel) {
 	/*
 	 * Segment nuclei using stardist, remove background and grow regions.
 	 * Segmented regions are added to the ROI manager
 	 * 
 	 */
 	print("Segmentation [ mem:" + round(IJ.currentMemory()/1e6)+"MB]");
 	name = getTitle();
 	selectImage(id);
 	run("Select None");
 	print(name); 	

 	print("Creating  a mask [ mem:" + round(IJ.currentMemory()/1e6)+"MB]"); 	 
 	mask = createBackgroundMask(id, 10); 	 	
 	 	 	
 	print("Running StarDist [ mem:" + round(IJ.currentMemory()/1e6)+"MB]");
 	run("Select None");
 	selectImage(id);
 	run("Duplicate...", "title=dapi duplicate channels="+channel);
 	nuc = getImageID();
 	run("Square Root");
 	run("Command From Macro", "command=[de.csbdresden.stardist.StarDist2D], args=['input':'dapi','modelChoice':'DSB 2018 (from StarDist 2D paper)', 'normalizeInput':'true', 'percentileBottom':'5.0', 'percentileTop':'95.0', 'probThresh':'0.5', 'nmsThresh':'0.4', 'outputType':'Label Image', 'nTiles':'9', 'excludeBoundary':'0', 'roiPosition':'Automatic', 'verbose':'false', 'showCsbdeepProgress':'false', 'showProbAndDist':'false'], process=[false]");
	selectImage("Label Image"); label = getImageID(); 	
	run("Remap Labels");
	
 	print("Running Watershed [ mem:" + round(IJ.currentMemory()/1e6)+"MB]"); 
 	selectImage(label); 		
	run("Label Size Filtering", "operation=Greater_Than size=100");	
	rename("Input");
	input = getImageID();
	
 	run("Marker-controlled Watershed", "input=[Input] marker=[Label Image] mask=[Mask] compactness=0 calculate use");
 	selectImage("Input-watershed"); watershed = getImageID();
 	run("Remap Labels");
 	run("Collect Garbage");	
 	
 	print("Creating ROIs [ mem:" + round(IJ.currentMemory()/1e6)+"MB]");
 	run("Select None"); 
 	amax = getValue("Area");
 	rois = createROIfromLabels(100, amax); 	
 	
 	closeImageList(newArray(nuc,mask,label,input,watershed)); 	
 	run("Collect Garbage");
 	
 	return rois;
}

function createBackgroundMask(id, percentile) {
	/*
	 * Create a background mask
	 */
	selectImage(id);
	run("Duplicate...", "title=tmp duplicate");		
	equalizeChannels(getImageID());	
	run("8-bit");
	id1 = getImageID();	
	run("Z Project...", "projection=[Average Intensity]");		
	run("Median...", "radius=5");	
	run("32-bit");
	t = calculatePercentile(percentile);
	setThreshold(t, getValue("Max"), "raw");
	run("Convert to Mask");
	id2 = getImageID();
	rename("Mask");	
	run("Grays");
	roiManager("select", 0);
	if (!isSelectionAllImage()) {
		run("Make Inverse");
		setColor(0);
		fill();
	}	
	selectImage(id1);close();	
	return id2;
}

function isSelectionAllImage() {
	/*return true if the selection is the all image or there is no selection*/	
	getPixelSize(unit, pixelWidth, pixelHeight);
	getDimensions(width, height, channels, slices, frames);
	if (getValue("Area") == (width * height * pixelWidth * pixelHeight)) {
		return true;
	} else {
		return false;
	}	
}

function calculatePercentile(percentage) {
	/* Compute a threshold based on a percentile of the image intensity */	
	if (bitDepth()==8) {
		nbins= 255;
		getHistogram(values, counts, nbins);
	} else {		
		minbins = getValue("Min");
		maxbins = getValue("Max");
		nbins = maxOf(100, maxbins - minbins);
		getHistogram(values, counts, nbins, minbins, maxbins);
	}
	N = 0;
	for (i = 0;	i < counts.length; i++) {
		N += counts[i];
	}
	n = 0;
	for (i = 0; i < counts.length && n < percentage / 100 * N; i++) {
		n += counts[i]; 
	}
	return values[minOf(i, nbins-1)];
}

function calculateThreshold() {
	/* Compute a threshold based on a percentile of the image intensity */
	
	return getValue("Median") + getValue("StdDev");		
}

function createROIfromLabels(area_min, area_max) {
	/*
	 * Create ROI from labels
	 */
	nlabels = getValue("Max");
 	print("Number of regions " + nlabels); 	
 	rois = newArray(nlabels);
 	n = 0;
 	for (i = 1; i <= nlabels; i++) {
 		run("Select None");
 		setThreshold(i,i);
 		run("Create Selection");
 		area = getValue("Area");
 		if (area > area_min && area < area_max) { 			
 			roiManager("add");
 			roiManager("select", roiManager("count")-1);
 			roiManager("remove slice info");
 			roiManager("rename", "ROI-"+n+1); 			
 			rois[n] = roiManager("count")-1;
 			n++;
 		}
 	}
 	return Array.trim(rois, n);
}

function measureInROI(tbl, id, rois, channel_names) {
	/*
	 * Measurement for each ROI over all chanels.
	 * - Mean
	 * - Area
	 * - Colocalization
	 */
	 print("Measure in ROIs");	 
	 Table.create(tbl1);
	 Stack.getDimensions(width, height, channels, slices, frames);
	 	 
	 for (i = 0; i < rois.length; i++) {
	 	roiManager("select", rois[i]);
	 	Table.set("Index", i, Roi.getName);
	 	
	 	// Area
	 	Table.set("Area", i, getValue("Area"));
	 	
	 	// Mean in each channel
	 	for (c1 = 1; c1 <= channels; c1++) {
	 		Stack.setChannel(c1);
	 		Table.set("Mean ["+channel_names[c1-1]+"]", i, getValue("Mean"));	 	 		
	 	}
	 	
	 	// colocalization	 	
	 	for (c1 = 1; c1 <= channels; c1++) {
			Stack.setChannel(c1);
			I1 = getIntensities();
			for (c2 = c1+1; c2 <= channels; c2++) {
				Stack.setChannel(c2);
				I2 = getIntensities();
				pcc = pearson(I1, I2);				
				Table.set("PCC ["+channel_names[c1-1]+":"+channel_names[c2-1]+"]", i, pcc[0]);				
			}
	 	}
	 }
	 Table.update;
	 run("Collect Garbage");	
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

function computeRobustThresholdOnArray(values, alpha) {
	tmp = Array.copy(values);
	Array.sort(tmp);
	m = tmp[round(tmp.length/2)];
	for (i= 0 ; i < tmp.length; i++) {
		tmp[i] = abs(tmp[i] - m);
	}
	Array.sort(tmp);
	return m + alpha * tmp[round(tmp.length/2)];	
}

function measureROIdistanceToLabels(tbl, id, rois) {
	selectImage(id);
	run("Options...", "iterations=1 count=1 edm=32-bit do=Nothing");	
	run("Distance Map", "stack");	
	Stack.getDimensions(width, height, nchannels, nslices, nframes);
	for (c = 1; c <= nchannels; c++) {
		for (i = 0; i < rois.length; i++) {
			roiManager("select", i);
			Stack.setChannel(c);
			Table.set("Distance to [" + channel_names[c-1]+"]", i, getValue("Min"));
		}
	}
	Table.update;
	close();
	run("Collect Garbage");	
}

function recordPositiveROIs(tbl, id, rois) {
	/*
	 * Fill the table with a flag telling if the cell is positive in this channel
	 * Parameters:
	 * tbl: name of the table
	 * id: index of the map
	 * rois: list of ROIs
	 * 
	 * Returns: nothing
	 */
	selectImage(id);
	Stack.getDimensions(width, height, nchannels, nslices, nframes);
	for (c = 1; c <= nchannels; c++) {
		for (i = 0; i < rois.length; i++) {
			roiManager("select", i);
			Stack.setChannel(c);
			v = getValue("Mean");
			if (v > 128) {				
				Table.set("Positive [" + channel_names[c-1]  + "]", i, 1);
			} else {
				Table.set("Positive [" + channel_names[c-1] + "]", i, 0);
			}
		}
	}
	Table.update;
}

function mapPostiveROIs(id, rois) {
	/*
	 * Create a map of positive ROI using a background sutraction and Otsu threshold 
	 * on the image of the median.
	 * 
	 * id : image ID
	 * rois : array of ROI id
	 * 
	 * Returns: the id of the map Image
	 */
	
	run("Select None");
	run("Duplicate...", "title=positive duplicate");
	Stack.setDisplayMode("grayscale");		
	run("Subtract Background...", "rolling=50 stack");		
	id1 = getImageID();		
	Stack.getDimensions(width, height, nchannels, nslices, nframes);
	for (c = 1; c <= nchannels; c++) {
		for (i = 0; i < rois.length; i++) {
			roiManager("select", rois[i]);
			run("Enlarge...", "enlarge=-1");
			Stack.setChannel(c);
			v = getValue("Median");
			setColor(v);
			fill();
		}
		roiManager("select", rois);	
		roiManager("Combine");
		run("Make Inverse");
		run("Enlarge...", "enlarge=-1");
		setColor(0); fill();
		setMinAndMax(0, 255);
	}
	run("Convert to Mask", "method=Otsu background=Dark calculate");	
	return id1;
}

function summarizeTable(src, dst, filename, channel_names) {
	/*
	 * Summarize the result in a table with a line per image
	 */
	if (!isOpen(dst)) {Table.create(dst);}
	selectWindow(dst);
	row = Table.size;
	Table.set("File", row, filename);
	// total number of ROI
	selectWindow(src);		
	n = Table.size;
	selectWindow(dst);
	Table.set("Number of ROIs", row, n);
	
	// count positive
	for (c = 0; c < channel_names.length; c++) {
		selectWindow(src);			
		x = Table.getColumn("Positive [" + channel_names[c]  + "]");
		n = 0;
		for (i = 0; i < x.length; i++) {
			if (x[i]==1) {
				n++;
			}
		}
		selectWindow(dst);
		Table.set("Positive [" + channel_names[c]  + "]", row, n);
	}
	
	// count double positive
	for (c1 = 0; c1 < channel_names.length; c1++) {
		for (c2 = 0; c2 < channel_names.length; c2++) if (c1 != c2) {
			selectWindow(src);					
			x = Table.getColumn("Positive [" + channel_names[c1]  + "]");
			y = Table.getColumn("Positive [" + channel_names[c2]  + "]");
			n = 0;
			for (i = 0; i < x.length; i++) {
				if (x[i]==1 && y[i]==1) {
					n++;
				}
			}
			selectWindow(dst);		
			Table.set("Double Positive [" + channel_names[c1]  + ":" + channel_names[c2] + "]", row, n);
		}
	}
	
	// average distance to other labels if positive
	for (c1 = 0; c1 < channel_names.length; c1++) {
		for (c2 = 0; c2 < channel_names.length; c2++) if (c1 != c2) {
			selectWindow(src);					
			x = Table.getColumn("Positive [" + channel_names[c1]  + "]");
			y = Table.getColumn("Distance to [" + channel_names[c2]  + "]");
			n = 0;
			d = 0;
			for (i = 0; i < x.length; i++) {
				if (x[i]==1) {
					n++;
					d += y[i];
				}
			}
			selectWindow(dst);					
			Table.set("Average distance [" + channel_names[c1]  + ":" + channel_names[c2] + "]", row, d/n);
		}
	}
	Table.update;
}

function cropActiveSelection() {	
	print("ROI size " + Roi.size);
	if (Roi.size == 0) {
		run("Select All");
		roiManager("add");
	} else {
		run("Crop");
		roiManager("add");
	}	
	roiManager("select", 0);
	roiManager("rename", "Selection");
}

function parseCSVString(csv) {
	str = split(csv,",");
	values = newArray(str.length);
	for (i = 0 ; i < str.length; i++) {
		values[i] = String.trim(str[i]);
	}
	Stack.getDimensions(width, height, channels, slices, frames);
	for (i = str.length; i < channels; i++) {
		values[i] = "ch" + i;
	}
	return values;
}

function equalizeChannels(id) { 	
 	/* Equalize all channels */
 	selectImage(id);
 	 run("Select None");
 	for (c = 1; c <= nSlices; c++) {
 		Stack.setChannel(c);
 		run("Enhance Contrast", "saturated=0.15");
 	}
 }
 
function closeWindow(name) {
 	if (isOpen(name)) { 		
 		selectWindow(name); 
 		run("Close");
 	}
}

function closeImageList(list) {
	for (i = 0 ;i < list.length; i++) {
		selectImage(list[i]);
		close();
	}
}

function addROIsToOverlay(id, rois) {
	selectImage(id);
	n = roiManager("count");
	for (i = 0; i < n; i++) {
		roiManager("select", i);
		run("Add Selection...");	
		Overlay.setStrokeColor("white");
	}		
}
