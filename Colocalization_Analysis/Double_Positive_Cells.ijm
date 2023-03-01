#@File    (label="Input file",description="Small tiff with ROI in it") filename
#@Integer (label="Reference Channel",description="Channel to segment",value=2) ch1
#@Float   (label="Sensitivity",description="Control the threshold to detect cells (dafulat:3)",value=3) sensitivity1
#@Float   (label="Cell size [um]", description="Minimum diameter of detected cells",value=10) minsz1
#@Integer (label="Measure Channel",description="Channel for intensity measure",value=1) ch2
#@Float   (label="Sensitivity",description="Control the threshold to detect double positive cells (default:3)",value=3) sensitivity2
#@OUTPUT ratio

/*
 * Segment and count double positive cells
 * 
 * Open and run the macro from the script editor
 * 
 * Tested on 1.53v
 * 
 * Jerome Boulanger 2022 for Elena
 */
 
setBatchMode("hide");
run("Close All");
if (isOpen("ROI Manager")) {selectWindow("ROI Manager"); run("Close");}

loadImage(filename);

print("Segmenting image");
process_image(ch1, ch2, sensitivity1, minsz1);


print("Measuring intensities in ROIs");
ratio = measure(filename, ch1, ch2, sensitivity2);

print("Saving visualization in file\n"+ replace(filename,".small.tif",".process.tif"));
run("Select None");
Stack.setChannel(ch1);
setMinAndMax(0, 6 * getValue("StdDev"));
Stack.setChannel(ch2);
setMinAndMax(0, 6 * getValue("StdDev"));
Stack.setDisplayMode("Composite");
saveAs("TIFF", replace(filename,".small.tif",".processed.tif"));

run("Close All");
if (isOpen("ROI Manager")) {selectWindow("ROI Manager"); run("Close");}
setBatchMode("exit and display");
print("Done\n\n");

function loadImage(filename) {
	/*
	 * Load small image and its ROI
	 * Transfer and scale the ROI to the full res image 	 
	 */	 
	print("Loading low res image " +filename);
	open(filename);	
	id1 = getImageID();
	w1 = getWidth();
	h1 = getHeight();
	if(!isThereAnActiveSelection()) {exit("No ROI in low res image");}
	roiManager("add");
	close();		

	full = replace(filename, ".small.tif", ".full.tif");	
	print("Loading full res image " + full);
	if (!File.exists(full)) {exit("Full resolution image not found")};
	open(full);	
	id2 = getImageID();
	w2 = getWidth();
	h2 = getHeight();
	print("Scaling selection");
	run("Restore Selection");	
	roiManager("select", 0);
	run("Scale... ", "x="+w2/w1+" y="+h2/h1);
	run("Interpolate", "interval=5 smooth");
	run("Duplicate...", "duplicate");
	roiManager("rename", "Selection");
	roiManager("update");
	//roiManager("delete");
	
}

function isThereAnActiveSelection() {
	/*
	 * Test if there is an active selection
	 */
	getPixelSize(unit, pixelWidth, pixelHeight);
	a = getValue("Area") / pixelWidth / pixelHeight;
	b = getWidth() * getHeight();	
	if (a == b) {
		return false;
	} else {
		return true;
	}
}

function process_image(ch1, ch2, sensitivity1, minsz1) {
	process_ref_channel(ch1, sensitivity1, minsz1);
	preprocess(ch1);
	preprocess(ch2);	
}

function process_ref_channel(ch, sensitivity, minsz) {
	/*
	 * Segment objects
	 * 
	 * ch: channel to segment
	 * sensivity : control teh segmentation threshold
	 * minsz : minimum diameter of the objects
	 * 
	 */
	print("Processing reference channel " + ch);
	run("Duplicate...", "duplicate channels="+ch);	
	run("Median...", "radius=2");
	run("Subtract Background...", "rolling=50");	
	setThreshold(getValue("Mean") + sensitivity * getValue("StdDev"), getValue("Max"));
	run("Convert to Mask");
	run("Watershed");
	run("Minimum...", "radius=5");
	run("Maximum...", "radius=5");	
	roiManager("select", 0);
	run("Analyze Particles...", "size="+(PI*minsz*minsz/4)+"-Infinity add");
	close();
}

function preprocess(ch) {
	/*
	 * Preprocess the channel with background subtraction
	 * 
	 */
	print("Pre-processing channel " + ch);
	run("Select None");
	Stack.setDisplayMode("color");
	Stack.setChannel(ch);
	run("Subtract Background...", "rolling=50 slice");
	run("Subtract...", "value="+(getValue("Mean")+getValue("StdDev"))+" slice");
}

function measure(filename, ch1, ch2, sensitivity) {	
	/*
	 * Measure the intenisty and count the number of object where intensity
	 * in channel ch2 is above a threhosld, report results in two tables
	 * 
	 */
	Stack.setChannel(ch2);
	threshold = getValue("Mean") + sensitivity * getValue("StdDev");
	n = roiManager("count");
	if (!isOpen("rois.csv")) { Table.create("rois.csv");}
	selectWindow("rois.csv");
	row = Table.size;		
	ndoublepositive = 0;
	ncells = n-1;
	print("table size " + row);
	for (i = 0; i < n; i++) {		
		roiManager("select", i);
		if (!matches(Roi.getName, "Selection")) {
			Table.set("Filename", row, filename);
			Table.set("Cell", row, i+1);		
			Table.set("Area [um^2]", row, getValue("Area"));
			Stack.setChannel(ch1);		
			Table.set("Mean Channel " + ch1, row, getValue("Mean"));
			Stack.setChannel(ch2);			
			Table.set("Mean Channel " + ch2, row, getValue("Mean"));		
			if (getValue("Mean") > threshold) {
				ndoublepositive++;
				Overlay.addSelection("yellow");	
				Overlay.setPosition(ch2, 0, 0);				
			} else{
				Overlay.addSelection("green");	
				Overlay.setPosition(ch1, 0, 0);
			}
			row++;
		} else {
			Overlay.addSelection("white");
			Overlay.setPosition(0, 0, 0);		
		}
		//Overlay.setPosition(0, 0, 0);		
		Overlay.show();
	}	
	
	setFont("Arial", 50, "scale");
	
	setColor("green");
	Overlay.drawString("Number of positive cells in channel #"+ch1+" : " + (ncells - 1), 5, 55);
	
	setColor("yellow");
	Overlay.drawString("Number of double positive double cells in channel #"+ch1+" and channel #"+ch2+" : " + ndoublepositive, 5, 105);
	
	run("Scale Bar...", "width=500 height=10 thickness=10 font=50 color=White background=None location=[Lower Right] horizontal bold overlay");
	
	Table.update();
	if (!isOpen("image.csv")) { Table.create("image.csv");}
	selectWindow("image.csv");
	row = Table.size;
	Table.set("Filename", row, filename);
	Table.set("Number of cells in channel " + ch1, row, ncells);
	Table.set("Number of cells with signal in channel " + ch2, row, ndoublepositive);
	Table.set("Ratio", row, ndoublepositive / ncells);
	return ndoublepositive / ncells;
}
