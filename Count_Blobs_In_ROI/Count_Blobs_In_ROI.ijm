// @Integer(label="Channel", value=1) channel
// @Integer(label="Spot Size [nm]", value=300) spot_size_nm
// @Integer(label="Specificity", value=2) specificity

/*
 * Count the number of 2D blobs in ROIs
 * 
 * Apprends results in a table counts.csv 
 * 
 * Jerome for Marco 2021
 */
 


setBatchMode(true);
run("Select None");

// list the ROI indices in an array
rois = Array.getSequence(roiManager("count"));

// Segment blobs and get a indices of the ROI
blobs = segmentBlobs(channel,spot_size_nm/1000.,specificity);

// counts the number of blobs with in each ROI
counts = countsBlobsInRoi(rois, blobs);

// append the results to the table
appendResults("Counts.csv", counts, channel);

// add overlay with specific channel
addOverlays(blobs,channel);

// remove blobs from the ROI manager
cleanUpROIManager(blobs);

run("Select None");
print(channel);
Stack.setChannel(channel);
setBatchMode(false);

function addOverlays(blobs,channel) {
	// add overlay
	Overlay.remove();
	for (i = 0; i < blobs.length; i++) {
		roiManager("select", blobs[i]);		
		Overlay.addSelection("white");		
		Overlay.setPosition(channel,1,1);		
	}
	Overlay.show();
}

function cleanUpROIManager(blobs) {	
	for (i = blobs.length-1; i >= 0; i--) {
		roiManager("select", blobs[i]);
		roiManager("delete");
	}
}

function appendResults(tblname, counts, channel) {
	if (!isOpen(tblname)) {
		Table.create(tblname);
	} else {
		selectWindow(tblname);
	}
	row = Table.size;
	for (i = 0; i < counts.length; i++) {
		Table.set("Image", row, getTitle());
		Table.set("Channel", row, channel);
		Table.set("ROI",row, i+1);
		Table.set("Count",row, counts[i]);
		
		row++;
	}
}

function countsBlobsInRoi(rois,blobs) {
	// count the number of blobs for each ROI
	counts = newArray(rois.length);
	for (i = 0; i < rois.length; i++) {		
		for (j = 0; j < blobs.length; j++) {
			if (hasRoi(rois[i], blobs[j]) == true) {
				counts[i] = counts[i] + 1;
			}
		}
	}
	return counts;
}

function segmentBlobs(channel,size_um,specificity) {
	// segment blobs in an image		
		
	run("Duplicate...", "title=id1 duplicate channels="+channel);	
	getPixelSize(unit, pixelWidth, pixelHeight);			
	run("32-bit");
	id1 = getImageID;	
	run("Gaussian Blur...", "sigma=" + size_um / pixelWidth);	
	run("Duplicate...", "title=id2 duplicate");
	id2 = getImageID;
	run("Gaussian Blur...", "sigma="+2 * size_um / pixelWidth);	
	imageCalculator("Subtract 32-bit stack",id1,id2);
	selectImage(id2);close();
	
	threshold(-specificity);
	
	run("Watershed");
	n0 = roiManager("count");
	min_size =  PI * size_um * size_um;
	run("Analyze Particles...", "size="+min_size+"-Infinity  add");
	n1 = roiManager("count");
	rois = newArray(n1-n0);
	for (i=0;i<rois.length;i++){
		rois[i] = n0 + i;
		roiManager("select", rois[i]);
		Roi.setPosition(channel, 1, 1);
		roiManager("update");
	}
	close();
	roiManager("show none");
	return rois;
}


function hasRoi(r1, r2) {
	// return true if the R1 contrains R2
	ret = false;
	roiManager("select", r2);
	Roi.getCoordinates(x, y);
	roiManager("select", r1);	
	for (i = 0; i < x.length; i++) {
		if (Roi.contains(x[i], y[i])) {
			ret = true;
		}
	}
	return ret;
}

function threshold(logpfa) {
	 // Define a threshold based on a probability of false alarm (<0.5)
	 st = laplaceparam();
	 pfa = pow(10,logpfa);	 
	 if (pfa > 0.5) {
	 	t = st[0] + st[1]*log(2*pfa);
	 } else {
	 	t = st[0] - st[1]*log(2*pfa);
	 }
	 setThreshold(t, st[2]);
	 run("Convert to Mask");	 
}

function laplaceparam() {
	/* compute laplace distribution parameters 
	 *  Return center, scale and maximum value
	 */
	getHistogram(values, counts, 100);	
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
	return newArray(m, b/n, values[values.length-1]);
}




