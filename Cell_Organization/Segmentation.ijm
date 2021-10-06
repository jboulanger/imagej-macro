run("Remove Overlay");

segment(getImageID, 4, false, 20, 1, 10, true);
segment(getImageID, 1, true,   5, 1,  1, false);

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
	cmap = newArray("#AA00AA","#55AA55","#FF9999","#FFFF55");
	for (i = n0; i < n1; i++) {
		roiManager("select", i);
		Roi.setStrokeColor(cmap[ch-1]);
		Roi.setStrokeWidth(2);
	}
	selectImage(id1); close();
	print("   - "+(n1-n0+1)+" ROI added");
	selectImage(id);
	roiManager("show all without labels");
	return a;
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
	med = computeArrayQuantile(A,0.5);	
	mad = computeMeanAbsoluteDeviation(A,med);		
	t = med + logpfa * log(10) * mad;
	print("  - med = " + med  + " mad = " + mad + " t = " + t);
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
