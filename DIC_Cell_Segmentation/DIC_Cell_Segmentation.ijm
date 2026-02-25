#@File (label="Input file") filepath
#@String (label="Measurements",value="AR,Circ.,AR,Round",description="comma separated list of measurements.") metrics_str
#@Float(label="Threshold factor",value=0.1) threshold_factor
#@Float(label="Min Area",value=100) min_area

/*
 * Segment DIC images from incucyte
 * 
 * 
 * Jerome Boulanger for David Paul - 2026
 */

//filepath =  "/home/jeromeb/work/userdata/McMahon_group/David Paul/DIC/raw_D2_3_04d00h00m.tif";

function solvePoissonFFT(regularization) {
	// Solve Laplacian x = Image	
	id = getImageID();	
	run("Duplicate...", "title=flt "); flt = getImageID();	
	run("Macro...", "code=[v=-1/(4+"+pow(10,-regularization)+"-2*cos(2*(x-w/2)*3.14159/w)-2*cos(2*(y-h/2)*3.14159/h))]");	
	selectImage(id);
	run("Custom Filter...", "filter=flt process");
	rename("Reconstructed"); rec= getImageID();
	selectImage(flt); close();
	selectImage(rec);		
}

function segmentDIC(id, threshold_factor, min_area) {
	/*Segment the DIC image*/
	selectImage(id);
	run("Duplicate..."," ");
	run("32-bit");
	run("Gaussian Blur...", "sigma=1");
	solvePoissonFFT(1e-3);
	//run("Median...", "radius=1");
	threshold = getValue("Median") + threshold_factor*getValue("StdDev");
	setThreshold(threshold, getValue("Max"));
	run("Analyze Particles...", "size="+min_area+"-Infinity add");
	//close();
	selectImage(id);
	roiManager("Show All without labels");
}

function measureROI(id, name, metrics, tbl) {
	/*Measure rois*/
	selectWindow(tbl);
	nrows = Table.size;
	selectImage(id);
	n = roiManager("count");
	for (index=0;index<n;index++) {
		roiManager("select", index);
		Table.set("name", nrows+index, name);
		Table.set("ROI", nrows+index, Roi.getName);
		for (midx = 0 ; midx < metrics.length; midx++) {
			value = getValue(metrics[midx]);
			Table.set(metrics[midx], nrows+index, value);
		}
	}
}

function closeWindow(title) {
	/* close the window if opened */
	if (isOpen(title)) {
		selectWindow(title);
		run("Close");
	}
}

function parseCSVString(csv, default) {	
	if (matches(csv, "auto")) {
		values = default;
	} else {
		str = split(csv,",");
		values = newArray(str.length);
		for (i = 0 ; i < str.length; i++) {
			values[i] = String.trim(str[i]);
		}
	}
	return values;
}

// start of the script
tbl = "DIC_Shape_results";

// get the metrics
default_metrics = newArray("AR","Circ.","Area","Round");
metrics = parseCSVString(metrics_str, default_metrics);

// close all and open the picture
run("Close All");
closeWindow("ROI Manager");
name = File.getNameWithoutExtension(filepath);
open(filepath);
id = getImageID();

// create a table if needed
if (!isOpen(tbl)) {
	Table.create(tbl);
}

// segment
segmentDIC(id, threshold_factor, min_area);

// measure
measureROI(id, name, metrics, tbl);


