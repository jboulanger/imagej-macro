//@Integer(label="Channel 1",value=1) ch1
//@Integer(label="Channel 2",value=2) ch2
//@String(label="Threshold",choices={"Otsu","Huang","Max Entropy"}) auto_threshold

/*
 * Basic overlap test
 * 
 * Input 2D hyperstacks
 * 
 * Compute Manders coefficients as
 *  M1 = sum_{ch2>threshold} ch1 / sum ch1
 *  M2 = sum_{ch1>threshold} ch2 / sum ch2
 *  
 * Compute overlap ratios as
 *  R1 = sum_{ch1>threshold & ch2>threshold} ch1 / sum_{ch1>threshold } ch1
 *  R2 = sum_{ch1>threshold & ch2>threshold} ch2 / sum_{ch2>threshold } ch2
 *  
 *  
 * Usage:
 *  1. Open the macro and the image
 *  2. Press run
 *  
 *  
 * Jerome Boulanger for Antoine Dessus 2025
 */


function createTest() {
	n = 512;
	
	newImage("Test", "8-bit black composite", n, n, 2, 1,1);
	
	Stack.setChannel(1);
	d = 0.1*(1 + random);
	x = n * (random - d);
	y = n * (random - d);
	makeOval(x, y, n*d, n*d);
	run("Make Band...", "band=5");
	val1 = 128*(1+random);
	setColor(val1);
	fill();
	
	Stack.setChannel(2);
	delta = 10;
	makeOval(x+delta*random, y+delta*random, n*d, n*d);
	run("Make Band...", "band=5");
	val2 = 128*(1+random);
	setColor(val2);
	fill();
	
	run("Select None");
	run("Add Noise", "stack");

	print(val1, val2);
	return getImageID();
}

function createROI(id, ch, auto_threshold) {
	/*
	 * Create ROIs thresholding the channels
	 * 
	 */
	print("Create ROI for channel "+ ch + " with threshold "+ auto_threshold);
	// Select the image
	selectImage(id);
	run("Select None");
	// Duplicate the channel
	run("Duplicate...", "duplicate channels=" + ch);
	// Get the ID of the channel image
	tmp = getImageID();
	// Apply a small median filter
	run("Median...", "radius=2");
	// Define the threshold
	setAutoThreshold(auto_threshold  + " dark");
	// Convert the thresholded area to a ROI
	run("Create Selection");
	// Add the ROI to the ROI manager
	roiManager("add");
	// Get the index of the added ROI
	roi = roiManager("count") - 1;
	// Close the temporary image
	selectImage(tmp);
	close();
	return roi;
}

function computeMandersCoefficient(id, ch1, ch2, roi1, roi2) {
	/*
	 * Compute manders coefficients M1 and M2
	 * 
	 */
	
	// Select the image
	selectImage(id);
	run("Select None");
	
	// Measure the intensity of ch1 in roi2
	roiManager("select", roi2);
	Stack.setChannel(ch1);
	sum_ch1_roi2 = getValue("RawIntDen");
	
	// Measure the intensity of ch2 in roi1
	roiManager("select", roi1);
	Stack.setChannel(ch2);
	sum_ch2_roi1 = getValue("RawIntDen");
	
	// Measure the intensity of ch1 
	run("Select None");
	Stack.setChannel(ch1);
	sum_ch1 = getValue("RawIntDen");
	
	// Measure the intensity of ch1 
	run("Select None");
	Stack.setChannel(ch2);
	sum_ch2 = getValue("RawIntDen");
	
	return newArray(sum_ch1_roi2 / sum_ch1, sum_ch2_roi1 / sum_ch2);
}


function computeOverlapCoefficient(id, ch1, ch2, roi1, roi2, overlap) {
	/*
	 * Compute overlap coefficients
	 * 
	 */

	// Select the image
	selectImage(id);

	// Measure the intensity of ch1 in roi1
	roiManager("select", roi1);
	Stack.setChannel(ch1);
	sum_ch1_roi1 = getValue("RawIntDen");
	
	// Measure the intensity of ch1 in overlap
	roiManager("select", overlap);
	Stack.setChannel(ch1);
	sum_ch1_overlap = getValue("RawIntDen");
	
	// Measure the intensity of ch2 in roi2
	roiManager("select", roi2);
	Stack.setChannel(ch2);
	sum_ch2_roi2 = getValue("RawIntDen");
	
	// Measure the intensity of ch2 in roi1
	roiManager("select", overlap);
	Stack.setChannel(ch2);
	sum_ch2_overlap = getValue("RawIntDen");

	return newArray(sum_ch1_overlap / sum_ch1_roi1, sum_ch2_overlap / sum_ch2_roi2);
}



if (nImages == 0) {
	id = createTest();
} else {
	id = getImageID();
}

// Close the ROI manager if opened
if (isOpen("ROI Manager")) {
	selectWindow("ROI Manager");
	run("Close");
}

run("Select None");
setBatchMode("hide");
// Create regions of interest in both channels
roi1 = createROI(id, ch1, auto_threshold);
roi2 = createROI(id, ch2, auto_threshold);
roiManager("select", newArray(roi1, roi2));
roiManager("and");
roiManager("add");
overlap = roiManager("count") - 1;

manders = computeMandersCoefficient(id,ch1,ch2,roi1,roi2);
ratios = computeOverlapCoefficient(id,ch1,ch2,roi1,roi2,overlap);

selectImage(id);
name = getTitle();
getPixelSize(unit, pixelWidth, pixelHeight);
getDimensions(width, height, channels, slices, frames);

// report results in a table
row = nResults;
setResult("Image name", row, name);
setResult("Width ["+unit+"]", row, width * pixelWidth);
setResult("Height ["+unit+"]", row, height * pixelHeight);
setResult("Pixel size ["+unit+"]", row, pixelWidth);
setResult("Channel 1", row, ch1);
setResult("Channel 2", row, ch2);
setResult("Threshold", row, auto_threshold);
setResult("M1", row, manders[0]);
setResult("M2", row, manders[1]);
setResult("Ratio 1", row, ratios[0]);
setResult("Ratio 2", row, ratios[1]);
setBatchMode("exit and display");