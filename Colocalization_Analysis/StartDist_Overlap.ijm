#@String(label="Channel Names", description="List names of the channels separated by coma (,)",value="Blue,Green,Red") channel_names_str
#@String(label="ROI colors", description="List colors for overlays separated by coma (,)",value="blue,green,red") roi_colors_str
#@Integer(label="NUmber of tiles", description="Number of tiles used in startdist segmentation",value=1) ntiles
/*
 * Count the number of cell body and their overlap
 * 
 * 
 * - Detect in each image nuclei-like strutures using stardist
 * - Comupute the number of overlapping ROIs 
 * 
 * Need a patch of the original stardist 
 * copy this file to the plugin folder in FIJI:
 * https://github.com/jboulanger/stardist-imagej/releases/download/0.3.0-noclear/StarDist_-0.3.0.jar
 * 
 * Jerome Boulanger for Elena Garcia 2022
 */

  
if (nImages==0) {
	createTestImage();
	channel_names_str = "Red,Green,Blue";
	roi_colors_str = "#FF5555,#55FF55,#55555FF";
}
start_time = getTime();
setBatchMode("hide");
channel_names = parseCSVString(channel_names_str);
roi_color_names = parseCSVString(roi_colors_str);
processImage(channel_names,roi_color_names);
setBatchMode("exit and display");
print("Finished in " + (getTime - start_time) / 1000 + " seconds.\n");

function parseCSVString(csv) {
	// split and remove trailing white spaces
	items = split(csv, ",");
	for (i = 0; i < items.length; i++) { 
		items[i] =  String.trim(items[i]);
	}
	return items;
}

function processImage(channel_names,roi_colors_names) {
	// process the image slice by slice	
	id = getImageID();
	Overlay.remove;
	Stack.getDimensions(width, height, channels, slices, frames);
	print("Image size " + width + "x" + height);
	if (!isOpen("Overlap")) {Table.create("Overlap");}	
	for (z = 1; z <= slices; z++) {			
		print("*** Processing slice " + z + "/" + slices + " ***");
		processSlice(id, z, channel_names, roi_colors_names);
	}		
}

function  processSlice(id, z, channel_names, roi_colors_names) {
	// process the z slice and report the results in a table
	closeRoiManager();	

	counts = newArray(channel_names.length); // number of ROI in each channel
	roi_subset = newArray();
	roi_channels = newArray(); // channel indices of each ROI
	for (c = 0; c < channel_names.length; c++) {
		print(" - segment channel " + channel_names[c]);
		roi = segmentPlane(id, z, c+1);
		roi_subset = Array.concat(roi_subset, roi);
		counts[c] = roi.length;
		colorRoi(roi, roi_colors_names[c], c+1, z);
		roi2Overlay(roi,z);
		Array.fill(roi, c);		
		roi_channels = Array.concat(roi_channels, roi);
		print("   count : " + counts[c]);
	}
	
	//print("ROI Indices");
	//Array.print(roi_subset);
	
	//print("ROI Channels");
	//Array.print(roi_channels);
	print("Computing overlap");
	selectImage(id);
	name =  getTitle();
	
	selectWindow("Overlap");
	n = Table.size;		
	Table.set("name",n,name);
	Table.set("slice",n,z);
	for (c = 0; c < channel_names.length; c++) {	
		Table.set("count in "+channel_names[c], n, counts[c]);
	}
	for (c1 = 0; c1 < channel_names.length; c1++) {		
		roi1 = findMatchInArray(roi_channels, c1);		
		roi1 = composeIndices(roi_subset, roi1);
		for (c2 = 0; c2 < channel_names.length; c2++) if (c1 != c2){
			roi2 = findMatchInArray(roi_channels, c2);
			roi2 = composeIndices(roi_subset, roi2);			
			res = overlap2(roi1, roi2);
			Table.set("overlap "+channel_names[c1]+":"+channel_names[c2], n, res);		
		}
	}	
	Table.update;
	closeRoiManager();
	run("Select None");
}
 
function composeIndices(subset, indices) {
	// return subset[indices]
	dst = newArray(indices.length);
	for (i = 0; i < indices.length; i++) {
		dst[i] = subset[indices[i]];
	}
	return dst;
}
 
function findMatchInArray(array, value) {
	// return an array with the indices of array matching value
	dst = newArray(array.length);
	k = 0;
	for (i = 0; i < array.length; i++) {
		if (abs(array[i]-value) < 0.01) {
			dst[k] = i;
			k++;
		}
	}
	return Array.trim(dst,k);	
}

function overlap(roi1,roi2) {	
	//print("ROI1");
	//Array.print(roi1);
	//print("ROI2");
	//Array.print(roi2);
	
	run("Select None");
	wxh = getValue("Width")*getValue("Height");	
	noverlap = 0;
	and = newArray(2);
	for (i = 0; i < roi1.length; i++) {
		for (j = 0; j < roi2.length; j++) if (roi1[i]!=roi2[j]){					
			and[0] = roi1[i];
			and[1] = roi2[j];			
			roiManager("select", and);			
			roiManager("and");				
			a = getValue("Area");			
			if ( a < wxh/2 && a > 1) {				
				noverlap++;
			}
		}
	}
	return noverlap;
}


function overlap2(roi1,roi2) {	
	//print("ROI1");
	//Array.print(roi1);
	//print("ROI2");
	//Array.print(roi2);	
	run("Select None");	
	noverlap = 0;	
	for (i = 0; i < roi1.length; i++) {
		roiManager("select", roi1[i]);
		Roi.getCoordinates(x, y);
		for (j = 0; j < roi2.length; j++) if (roi1[i]!=roi2[j]) {
			roiManager("select", roi2[j]);			
			for (k = 0; k < x.length; k++) {
				if (Roi.contains(x[k], y[k])) {
					noverlap++;
					break;				
				}
			}
		}
	}
	return noverlap;
}

function segmentPlaneLabel(id,z,c) {
	//use the label output from stardist plugin to get the ROIs
	selectImage(id);
	run("Select None");
	run("Duplicate...", "title=mask duplicate channels="+c+" slices="+z);
	id1 = getImageID();			
	run("Command From Macro", "command=[de.csbdresden.stardist.StarDist2D], args=['input':'mask', 'modelChoice':'Versatile (fluorescent nuclei)', 'normalizeInput':'true', 'percentileBottom':'0.0', 'percentileTop':'100.0', 'probThresh':'0.479071', 'nmsThresh':'0.3', 'outputType':'Label Image', 'nTiles':'"+ntiles+"', 'excludeBoundary':'2', 'roiPosition':'Automatic', 'verbose':'false', 'showCsbdeepProgress':'false', 'showProbAndDist':'false'], process=[false]");
	wait(5000);	
	//selectWindow("Label Image");	
	selectWindow("Label Image");
	id2 = getImageID();
	rois = label2ROI(100);	
	selectImage(id2); close();
	selectImage(id1); close();
	//setBatchMode(true);
	selectImage(id);
	return rois;	
}

function segmentPlane(id,z,c) {
	// output directly the ROI from stardist but startdist clears up the ROI manager (??)
	selectImage(id);
	run("Select None");
	run("Duplicate...", "title=mask duplicate channels="+c+" slices="+z);
	id1 = getImageID();	
	n0 = roiManager("count");
	run("Command From Macro", "command=[de.csbdresden.stardist.StarDist2D], args=['input':'mask', 'modelChoice':'Versatile (fluorescent nuclei)', 'normalizeInput':'true', 'percentileBottom':'0.0', 'percentileTop':'100.0', 'probThresh':'0.479071', 'nmsThresh':'0.3', 'outputType':'ROI Manager', 'nTiles':'"+ntiles+"', 'excludeBoundary':'2', 'roiPosition':'Automatic', 'verbose':'false', 'showCsbdeepProgress':'false', 'showProbAndDist':'false'], process=[false]");
	wait(500); // wait a bit
	n1 = roiManager("count");	
	
	// tidy up
	selectImage(id1); close();
	selectImage(id);
	
	// List the indices of the recently added ROI whose area > 100um
	rois = newArray(n1-n0);
	k = 0;
	for (i = 0; i < rois.length; i++) {
		roiManager("select", n0 + i);
		a = getValue("Area");		
		if (a > 20) {
			rois[k] = n0 + i;
			k++;
		}
	}
	rois = Array.trim(rois, k);
	return rois;	
}


function label2ROI(area_min) {
	run("Select None");
	wxh = getValue("Width")*getValue("Height");	
	getStatistics(area, mean, min, max, std, histogram);	
	count = 0;
	indices = newArray(max);
	print("Number of levels : " + max);
	for (i = 1; i <= max; i++) {		
		run("Select None");		
		setThreshold(i-0.1, i+0.1);
		run("Create Selection");
		a = getValue("Area");
		if ( a > area_min && a < wxh) {							
			roiManager("Add");
			indices[count] = roiManager("count") - 1;
			count++;
		}
	}
	print("counts : " + count);
	run("Select None");
	indices = Array.trim(indices, count);
	return indices;
}

function segmentPlaneold(z,c,scale,lambda) {	
	id = getImageID();
	run("Select None");
	run("Duplicate...", "duplicate channels="+c+" slices="+z);	
	scale=5;	lambda=5;
	run("Minimum...", "radius="+scale);
	run("Maximum...", "radius="+scale);	
	run("32-bit");
	run("Gaussian Blur...", "sigma="+scale+" stack");
	id1 = getImageID();
	run("Duplicate...", "use");
	run("Gaussian Blur...", "sigma="+3*scale+" stack");
	id2 = getImageID();
	imageCalculator("Subtract", id1, id2);
	selectImage(id2); close();
	selectImage(id1);
	getStatistics(area, mean, min, max, std, histogram);
	setThreshold(mean + lambda * std, max);
	run("Convert to Mask");
	run("Watershed");
	n0 = roiManager("count");
	run("Analyze Particles...", "size="+PI*scale*scale/4+"-Infinity add");	
	n1 = roiManager("count");	
	/*
	for (i = n1-1; i >= n0; i--) {
		roiManager("select", i);
		if (getValue("Area") < 10) {
			roiManager("delete");
		}
	}
	n1 = roiManager("count");
	*/
	ids = newArray(n1-n0);
	for (i = 0; i < ids.length; i++) {
		ids[i] = n0 + i;
	}
	selectImage(id1); close();	
	selectImage(id);
	roiManager("show none");
	return ids;
}


function closeRoiManager() {
	// close the roi manager
	while (roiManager("count")>0) {roiManager("select", roiManager("count")-1);roiManager("delete");}
	closeWindow("ROI Manager");
}

function closeWindow(name) {
	if (isOpen(name)) {
		selectWindow(name);
		run("Close");
	}
}

function colorRoi(ids,color,channel, slice) {
	for (i = 0; i < ids.length; i++) {
		roiManager("select", ids[i]);
		Roi.setStrokeColor(color);
		Roi.setPosition(channel, slice, 1);
		roiManager("update");
	}
}

function roi2Overlay(ids,z) {
	for (i = 0; i < ids.length; i++) {
		roiManager("select", ids[i]);
		run("Add Selection...");
		Overlay.setStrokeWidth(5);
		Overlay.show();
		Overlay.setPosition(1, z, 1);
	}
}


function createTestImage() {
	w = 500;
	h = 500;
	newImage("Test Image", "8-bit composite", w, h, 3 ,1,1);
	n = 3; // number of line and row 
	setColor(255);
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			x = (0.1+0.8*i / (n-1) - 0.25/n) * (w - 2/n);
			y = (0.1+0.8*j / (n-1) - 0.25/n) * (h - 2/n);
			makeOval(x, y, w/(2*n), h/(2*n));
			for (c = 1; c <= 3; c++) {
				if (random < 0.5) {				
					setSlice(c);
					fill();
				}
			}
		}
	}
	run("Select None");
}



