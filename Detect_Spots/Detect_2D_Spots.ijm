// @Integer(label="Channel", value=1) channel
// @Double(label="Spot size [px]", value=1.0, style=slider, min=0.0, max=3.0) spot_size
// @Double(label="Probablity of false alarm (log10)", value=1, style=slider, min=0, max=5) pfa
// @Boolean(label="Render", value=true) dorender
// @Double(label="Rendering spot size [px]", value=0.5, style=slider, min=0.0, max=20.0) rendering_spot_size
//
//
// Detection of spots in 3D stacks
//
// Spots are localized with a difference of Gaussian (DoG)
// The local maxima of the DoG and a threshold based on a false alarm rate are combined 
// to localize bright local maximas. Localized spots are then rendered as Gaussian spots
// and added as a second channel for visual validation. 
// 
//
// Usage: 
// 1. Draw roi where spots need to be detected and add them to the ROI manager (press key [t])
// 2. Run the macro
//
// Jerome Boulanger 2020 for Alejandro 
//

if (nImages==0) {
	generateTestImage();
}
	
img = getImageID; 
title = getTitle();
Stack.getDimensions(width, height, channels, slices, frames);

if (slices==1) {
	spot = detect2DSpots(channel, spot_size, pfa);
	count2DSpotsPerROI(spot, title);
	if (dorender) {
		overlaySpotImage(img, channel, spot, rendering_spot_size);
	}
} else {
	spot = detect3DSpots(channel, spot_size, pfa);
	count3DSpotsPerROI(title);	
	if (dorender) {
		renderSpotImage(img, channel, spot, rendering_spot_size);
	}
}
selectImage(spot);close();


function count3DSpotsPerROI(name) {
	getVoxelSize(w, h, d, unit);
	voxelVolume = w*h*d;
	for (i = 0; i < roiManager("count"); i++) {		
		roiManager("select", i);		
		Stack.getStatistics(voxelCount, mean, min, max, stdDev);		
		n = round(mean * voxelCount / max);	
		row = nResults;
		setResult("Name", row, name);	
		setResult("ROI", row, i+1);
		setResult("Number of Spots", row, n);
		setResult("ROI volume ["+unit+"]", row, voxelCount * voxelVolume);
		setResult("Spot density [Spot/"+unit+"]", row, n/(voxelCount*voxelVolume));
		setResult("Spot size [px]", row, spot_size);
		setResult("pfa", row, pfa);
	}
	updateResults();
}

function count2DSpotsPerROI(spot, name) {
	selectImage(spot);
	getPixelSize(unit, w, h);	
	pixelArea= w*h;
	for (i = 0; i < roiManager("count"); i++) {		
		roiManager("select", i);						
		n = getValue("RawIntDen") / 255;
		area = getValue("Area");
		row = nResults;
		setResult("Name", row, name);	
		setResult("ROI", row, i+1);
		setResult("Number of Spots", row, n);
		setResult("ROI volume ["+unit+"]", row, area);
		setResult("Spot density [Spot/"+unit+"]", row, n/area);
		setResult("Spot size [px]", row, spot_size);
		setResult("pfa", row, pfa);
	}
	updateResults();
}

function detect3DSpots(channel, size, pfa) {
	setBatchMode(true);
	run("Select None");	
	p = round(2*size+1);
	id0 = getImageID;
	// compute a difference of Gaussian
	run("Duplicate...", "title=id1 duplicate channels="+channel);	
	id1 = getImageID;	
	run("32-bit");
	
	run("Gaussian Blur 3D...", "x="+sigma1+" y="+sigma1+" z="+sigma1);
	run("Duplicate...", "title=id2 duplicate");
	id2 = getImageID;
	run("Gaussian Blur 3D...", "x="+sigma2+" y="+sigma2+" z="+sigma2);	
	imageCalculator("Subtract 32-bit stack",id1,id2);
	selectImage(id2);close();
	
	// local maxima
	selectImage(id1);
	run("Duplicate...", "title=id3 duplicate");
	id3 = getImageID;	
	run("Maximum 3D...", "x="+p+" y="+p+" z="+p);	
	imageCalculator("Subtract 32-bit stack",id3,id1);
	setThreshold(-0.1, 0.1);
	run("Make Binary", "method=Default background=Dark black");	
	run("Divide...", "value=255 stack");
	
	// threshold
	selectImage(id1);
	Stack.getStatistics(voxelCount, mean, min, max, stdDev);
	lambda = -2*norminv(pow(10,-pfa));
	//print("pfa="+pfa+" ->" +pow(10,-pfa) + "->"+lambda);
	threshold = mean + lambda * stdDev;
	setThreshold(threshold, max);	
	run("Make Binary", "method=Default background=Dark black");	
	run("Minimum 3D...", "x=1 y=1 z=1");
	run("Maximum 3D...", "x=1 y=1 z=1");
	run("Divide...", "value=255 stack");	
	
	// combine local max and threshold
	imageCalculator("Multiply stack",id1,id3);
	selectImage(id3);close();
	setThreshold(0.5, 1);
	run("Make Binary", "method=Default background=Dark black");		
	
	// remove edges	
	setSlice(1);
	run("Select All");
	setBackgroundColor(0, 0, 0);
	run("Cut");
	run("Select None");
	setSlice(nSlices);
	run("Select All");
	run("Cut");
	run("Select None");

	// count
	Stack.getStatistics(voxelCount, mean, min, max, stdDev);	
	print("N = " + mean * voxelCount / max);
	
	rename("Spot detection");
	setBatchMode(false);
	return id1;
}

function detect2DSpots(channel, size, pfa) {
	setBatchMode(true);
	run("Select None");	
	channel=2;
	sigma1 = size;
	sigma2 = 3 * size;
	p = round(2*size+1);
	id0 = getImageID;
	// compute a difference of Gaussian
	run("Duplicate...", "title=id1 duplicate channels="+channel);	
	id1 = getImageID;	
	run("32-bit");
	run("Square Root");
	run("Gaussian Blur...", "sigma="+sigma1);	
	run("Duplicate...", "title=id2 duplicate");
	id2 = getImageID;
	run("Gaussian Blur...", "sigma="+sigma2);	
	imageCalculator("Subtract 32-bit stack",id1,id2);
	selectImage(id2);close();
	
	// local maxima
	selectImage(id1);
	run("Gaussian Blur...", "sigma="+sigma1);
	run("Duplicate...", "title=id3 duplicate");	
	id3 = getImageID;	
	run("Maximum...", "radius="+p);
	imageCalculator("Subtract 32-bit stack",id3,id1);
	setThreshold(-0.00001, 0.00001);
	run("Make Binary", "method=Default background=Dark black");	
	run("Divide...", "value=255 stack");
	
	// threshold
	selectImage(id1);
	getStatistics(area, mean, min, max, stdDev);
	lambda = -2*norminv(pow(10,-pfa));
	print("pfa="+pfa+" ->" +pow(10,-pfa) + "->"+lambda);
	threshold = mean + lambda * stdDev;
	setThreshold(threshold, max);	
	run("Make Binary", "method=Default background=Dark black");	
	run("Minimum...", "radius=1");
	run("Maximum...", "radius=1");
	run("Divide...", "value=255 stack");	
	
	// combine local max and threshold
	imageCalculator("Multiply stack",id1,id3);
	selectImage(id3);close();
	setThreshold(0.5, 1);
	run("Make Binary", "method=Default background=Dark black");		
		
	// count
	n = getValue("RawIntDen") / 255;
	print("N = " + n);
	
	rename("Spot detection");
	setBatchMode(false);
	return id1;
}

function norminv(x) {
	return -sqrt(  log(1/(x*x)) - 0.9*log(log(1/(x*x))) - log(2*3.14159) );
}

function renderSpotImage(img, channel, spot, rendering_spot_size) {
	selectImage(spot);
	spot_title = getTitle;
	if (rendering_spot_size > 0.1) {
		run("32-bit");	
		run("Gaussian Blur 3D...", "x="+rendering_spot_size+" y="+rendering_spot_size+" z="+rendering_spot_size);
	}
	resetMinAndMax();
	run("16-bit");
	
	selectImage(img);
	
	run("Duplicate...", "title=id4 duplicate channels="+channel);	
	run("16-bit");
	run("Merge Channels...", "c1=id4 c2=["+spot_title+"] create");
	rename("Spot rendering");
	Stack.setChannel(2); 
	resetMinAndMax();	
}

function overlaySpotImage(img, channel, spot, rendering_spot_size) {
	selectImage(spot);	
	setThreshold(128, 255);
	run("Create Selection");
	Roi.getContainedPoints(xpoints, ypoints);
	selectImage(img);
	Overlay.remove;
	for (i = 0; i < xpoints.length; i++) {
		Overlay.drawEllipse(xpoints[i]-rendering_spot_size/2, ypoints[i]-rendering_spot_size/2, rendering_spot_size,rendering_spot_size);
		Overlay.setPosition(0, 0, 0);
	}
	Overlay.show;	
	resetMinAndMax();	
}

function generateTestImage() {
	w = 100;
	h = 100;
	d = 1;
	n = 10;
	newImage("Stack", "8-bit grayscale-mode", w, h, 1, d, 1);
	for (i = 0; i < n; i++) {
		x = 10 + round((w-10) * random);
		y = 10 + round((h-10) * random);
		z = 10 + Math.ceil((d-10) * random);
		//setSlice(z);
		setPixel(x,y,255);
	}
	run("Maximum 3D...", "x=2 y=2 z=2");
	run("Gaussian Blur 3D...", "x=0.5 y=0.5 z=0.5");
	run("Enhance Contrast...", "saturated=0");
	run("Add Noise", "stack");
	if (isOpen("ROI Manager")) {selectWindow("ROI Manager");run("Close");}
	makeRectangle(0, 0, w/2, h/2);
	roiManager("Add");
	makeRectangle(w/2, 0, w/2, h/2);
	roiManager("Add");
	makeRectangle(0, h/2, w/2, h/2);
	roiManager("Add");
	makeRectangle(w/2, h/2, w/2, h/2);
	roiManager("Add");
	run("Select None");
}
