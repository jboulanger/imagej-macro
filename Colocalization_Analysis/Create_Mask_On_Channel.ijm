// @File(label="Input",description="Filename or type 'image' to process opened image",value="image") path
// @Integer(label="Channels", value=1, description="index of the channel to process") channel

/*
 * Segment neuronal cells in a channel and add the mask as extra channel
 * 
 * This macro can be used to create a mask for the macro 
 * MultiChannel_Spot_3D_Colocalization.ijm
 * 
 * Requirements
 * ------------
 * Depends on FeatureJ, morpholibJ
 * 
 * Activate the IJPB-plugins and ImageScience update sites
 * 
 * Jerome Boulanger for Heidy Chen 2024
 */

function getDataSet() {
	/* Get the dataset given the path parameter */
	if (matches(path, ".*image")) {	
		print("Processing " + getTitle());
		return getImageID();
	} else {
		print("Processing " + path);
		run("Close All");		
		//open(path);		
		run("Bio-Formats Importer", "open=["+path+"] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
		return getImageID();
	}
}

function createMaskPlane(id) {
	/* Create a mask in place plane by plane */
	selectImage(id);
	Property.set("CompositeProjection", "null");
	Stack.setDisplayMode("grayscale");
	Stack.setChannel(channel);
	Stack.getDimensions(width, height, channels, slices, frames);	
	m = 0;
	s = 0;
	for (n = 1; n <= slices; n++) {
		Stack.setSlice(n);
		for (iter = 0; iter < 5; iter++) { run("Median...", "radius=2 slice");}
		run("Subtract Background...", "rolling=50 slice");
		
		m += getValue("Mean");
		s += getValue("StdDev");
	}
	m /= slices;
	s /= slices;
	threshold = m + 0.5*s;
	print(threshold);
	
	for (n = 1; n <= slices; n++) {
		Stack.setSlice(n);
		run("Macro...", "code=v=(v>="+threshold+") slice");		
		run("Minimum...", "radius=2 slice");
		run("Median...", "radius=2 slice");
		run("Maximum...", "radius=5 slice");
	}
	setMinAndMax(0, 5);
	Stack.setDisplayMode("composite");
	Stack.setSlice(round(slices/2));
}


function smooth() {
	/* Image enhancement for neurons */
	run("Duplicate...", "duplicate");
	run("32-bit");
	for (i = 0; i < 3; i++) {
		img1 = getImageID();
		run("Square Root", "stack");
		run("FeatureJ Hessian", "smallest smoothing=2");	
		img2 = getImageID();
		Stack.getStatistics(voxelCount, mean, min, max, stdDev);
		run("Macro...", "code=v=(-v*((-v)>0)) stack");
		imageCalculator("Multiply stack", img2, img1);
		run("Median 3D...", "x=2 y=2 z=2");
		selectImage(img1); close();
	}
}

function createMask(id) {
	/* Create a new image segmenting the channel in 3D */
	name = getTitle();
	Stack.getDimensions(width, height, channels, slices, frames);
	run("Split Channels");
	selectWindow("C"+channel+"-"+name);
	
	id0 = getImageID();
	
	smooth();

	// Threshold the smoothed image 
	Stack.getStatistics(voxelCount, mean, min, max, stdDev);	
	setThreshold(stdDev, max);
	setOption("BlackBackground", true);
	run("Convert to Mask", "black");	
	run("Maximum 3D...", "x=5 y=5 z=2");	
	id1 = getImageID();
	
	// Keep larger region of interest
	run("Connected Components Labeling", "connectivity=6 type=[16 bits]");
	id2 = getImageID();
	
	run("Label Size Filtering", "operation=Greater_Than size=30000");
	id3 = getImageID();	
	setThreshold(1, 65366);
	run("Convert to Mask", "black");
	run("16-bit");
	setMinAndMax(0, 1000);
	
	themask = getTitle();
	
	// Merge back the channels
	merge = "";
	for (c = 1; c <= channels; c++) {
		if (c == channel) {
			merge = merge + "c"+c+"=["+themask+"] ";
		} else {
			merge = merge + "c"+c+"=[C"+c+"-"+name+"] ";
		}
	}	
	
	run("Merge Channels...", merge + " create");
	dest = getImageID();	
	selectImage(id0); close();
	selectImage(id1); close();
	selectImage(id2); close();
	return dest;
}

function saveAndFinalize(id, path) {	
	/* Save results if an image was opened */
	if (!matches(path, ".*image")) {
		folder = File.getDirectory(path);
		fname = File.getNameWithoutExtension(path);		
		selectImage(id);
		path = folder + File.separator + fname + "-with-mask.tif";
		print("Saving result in file:");
		print(path);
		saveAs("TIFF", path);
	}
}

function main() {
	setBatchMode("hide");
	start_time = getTime();
	id = getDataSet();
	mask = createMask(id);
	saveAndFinalize(mask, path);
	setBatchMode("exit and display");
	end_time = getTime();
	print("Elapsed time " + (end_time - start_time)/1000 + " seconds.");
}

main();