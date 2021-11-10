#@Integer (label="ROI channel", value=1) roi_channel
#@Float(label="Rolling ball radius", value=200) ball_radius
#@Float(label="Specificity",description="decrease this value to", value=0.25) specificity
#@Float(label="min area [um^2]", value=50) min_area

/*
 * Segment image using smoothing and rolling ball and populate ROI Manger
 * 
 * Jerome Boulanger for Marco Laub 23021
 */
if (isOpen("ROI Manager")) {selectWindow("ROI Manager");run("Close");}
roiManager("count");
setBatchMode(true);
id0 = getImageID();
run("Select None");
// select a channel
run("Duplicate...", "duplicate channels="+roi_channel);
id = getImageID();
// smooth the image with a median filters
for (i = 0; i < 4; i++) {	
	run("Median...", "radius=2");
	run("Maximum...", "radius=2");
	run("Minimum...", "radius=2");
}
// substract the background
run("Subtract Background...", "rolling="+ball_radius);
// threshold
defineThreshold(-specificity);
run("Convert to Mask");
run("Watershed");
// close the roi manager if it is open

// add regions to the roi manager
run("Analyze Particles...", "size="+min_area+"-Infinity add");
filterOutROI(newArray("Circ.","Solidity"),newArray(0.1,0.1),newArray(1.0,1.0));
// close the processed image
selectImage(id); close();
selectImage(id0);
Stack.setChannel(roi_channel);
roiManager("show all without labels");
setBatchMode(false);

function defineThreshold(logpfa) {
	 // Define a threshold based on a probability of false alarm (<0.5)	 
	 getHistogram(values, counts, 1000); 
	 st = laplaceparam();
	 pfa = pow(10,logpfa);
	 if (pfa > 0.5) {
	 	t = st[0] + st[1]*log(2*pfa);
	 } else {
	 	t = st[0] - st[1]*log(2*pfa);
	 }
	 setThreshold(t, st[2]);
}

function laplaceparam() {
	// compute laplace distribution parameters 
	getHistogram(values, counts, 1000);
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
	return newArray(m,b/n,values[values.length-1]);
}

function filterOutROI(feat,min,max) {
	// filter out roi based on feature lists
	Array.print(feat);
	Array.print(min);
	Array.print(max);
	n = roiManager("count");
	for (i = n-1; i > 0; i--) {		
		roiManager("select", i);
		flag = false;
		for (j = 0; j < feat.length; j++) {		
			v = getValue(feat[j]);	
			if (v < min[j] || v > max[j]) {
				flag == true;				
			}
		}
		if (flag) {roiManager("delete");}
	}
}
