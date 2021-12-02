// @Integer(label="Channel", value=1) channel
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
blobs = segmentBlobs(2);

// counts the number of blobs with in each ROI
counts = countsBlobsInRoi(rois, blobs);

// append counts in a table
if (!isOpen("Counts.csv")) {
	Table.create("Counts.csv");
} else {
	selectWindow("Counts.csv");
}
row = Table.size;
for (i=0;i<rois.length;i++) {
	Table.set("Image", row, getTitle());
	Table.set("ROI",row, i+1);
	Table.set("Count",row, counts[i]);
	row++;
}

// cleanup the roi manager and add overlays
Overlay.remove();
for (i = blobs.length-1; i >= 0; i--) {
	roiManager("select", blobs[i]);
	Overlay.addSelection("red");
	
	roiManager("delete");
}
Overlay.show();
run("Select None");
setBatchMode(false);

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

function segmentBlobs(channel) {	
	// segment blobs in an image
	run("Duplicate...", "duplicate channels="+channel);
	run("Subtract Background...", "rolling=10");
	run("Median...","radius=1");
	getStatistics(area, mean, min, max, std, histogram);
	setThreshold(mean+2*std, max);
	run("Convert to Mask");
	run("Watershed");
	n0 = roiManager("count");
	run("Analyze Particles...", "size=5-Infinity pixel add");
	n1 = roiManager("count");
	rois = newArray(n1-n0);
	for (i=0;i<rois.length;i++){
		rois[i] = n0 + i;
	}
	close();
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




