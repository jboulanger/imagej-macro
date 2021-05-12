// @Integer(label="channel", value=1) channel
// @Double(label="spot size", value=1.0, style=slider, min=0.0, max=3.0) spot_size
// @Double(label="10^-{probablity false alarm}", value=1, style=slider, min=0, max=5) pfa
// @Boolean(label="render", value=true) dorender
// @Double(label="rendering spot size", value=0.5, style=slider, min=0.0, max=3.0) rendering_spot_size
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
	
img = getImageID; 
title = getTitle();
spot = detect3DSpots(channel, spot_size, pfa);
countSpotsPerROI(title);
if (dorender) {
	renderSpotImage(img, channel, spot,rendering_spot_size);
}

function countSpotsPerROI(name) {
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
	}
	updateResults();
}

function detect3DSpots(channel, size, pfa) {
	setBatchMode(true);
	run("Select None");
	sigma1 = size;
	sigma2 = 2 * size;
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
