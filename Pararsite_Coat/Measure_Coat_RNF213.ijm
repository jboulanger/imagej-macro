//@File(label="Input",description="Use image to run on current image",value="image") path
//@Boolean(label="Use shell") use_shell
//@Float(label="start time [min]",value=23) start_time
//@Float(label="first frame [frame]",value=23) first_frame
//@Boolean(label="Correct bleaching in fit") correct_bleaching

/*
	Measure RNF213 coat fluorescence overtime
		
	Segment the parasites in channel 1 and measure mean
	intensity on a shell around the parasite in channel 2.
	
	We assume that there is only a parasite centered in the field of
	view at the start of the sequence.
	
	Installation:
	- Install morpholibJ
	- Open the macro using the script editor in Fiji.
		
	Usage:
	
	- Open the image to process
	- Duplicate the region of interest with only 1 parasite
	- Run the macro
	- Save the table	
	
	
	Jerome Boulanger for Ana Crespillo Casado 2023
*/

function findNearestLabel(labels_frame, position) {
	/*
	 * Find nearest label and update position
	 * Return id of the image with the correct label
	 */
	selectImage(labels_frame);
	name = getTitle();
	run("Analyze Regions 3D", "centroid surface_area_method=[Crofton (13 dirs.)] euler_connectivity=6");
	selectWindow(name+"-morpho");
	
	xk = Table.getColumn("Centroid.X");
	yk = Table.getColumn("Centroid.Y");
	zk = Table.getColumn("Centroid.Z");
	labs = Table.getColumn("Label");
	print("Found " + labs.length +  " labels on frame " + frame);
	
	dmin = 1e6;
	kstar = 0;
	for (k = 0; k < xk.length; k++) {
		dx = position[0] - xk[k];
		dy = position[1] - yk[k];
		dz = position[2] - zk[k];
		d = sqrt(dx*dx+dy*dy+dz*dz);
		if (d < dmin) {
			dmin = d;
			kstar = k;
		}
	}
	
	position[0] = xk[kstar];
	position[1] = yk[kstar];
	position[2] = zk[kstar];
	
	selectWindow(name+"-morpho");
	run("Close");
	
	selectImage(labels_frame);
	run("Select Label(s)", "label(s)=" + labs[kstar]);
	selectWindow(name+"-keepLabels");
	labels_keep = getImageID();
	run("Remap Labels");
		
	return labels_keep;
}


function segmentAndTrackParasite(img) {
	/*
	 * Segment and track the parasite
	 * 
	 * Input:
	 * 	img : input image id with channel 1 to be segmented
	 * 	
	 * Output:
	 * 	labels: labels window id
	 */
	if (isOpen("labels")) {
		selectWindow("labels");
		close();
	}
	Stack.getDimensions(width, height, channels, slices, frames);
	Stack.getUnits(X, Y, Z, Time, Value);
	getVoxelSize(dx, dy, dz, unit);	
	run("Duplicate...", "title=binary duplicate channels=1");
	for (iter = 0; iter < 5; iter++) {
		run("Median 3D...", "x=2 y=2 z=1");
	}
	Stack.getStatistics(voxelCount, mean, min, max, stdDev);
	setThreshold(mean+stdDev, max+1);	
	run("Convert to Mask", "background=Dark black");	
	id0 = getImageID();
	
	position = newArray(width / 2 * dx, height / 2 * dy, slices / 2 * dz);
	
	for (frame = 1; frame <= frames; frame++) {
		
		print("Processing frame " + frame);
		selectImage(id0);
		run("Duplicate...", "title=test duplicate frames="+frame);
		id1 = getImageID();
		
		run("Connected Components Labeling", "connectivity=6 type=[8 bits]");
		selectWindow("test-lbl");
		id2 = getImageID();	
		
		keeplabel = findNearestLabel(id2, position);
		
		selectImage(id2); close();
							
		if (frame > 1) {	
			run("Concatenate...", "title=labels open image1=labels image2=test-lbl-keepLabels image3=[-- None --]");		
			selectWindow("labels");
			labels = getImageID();								
		} else {
			selectImage(keeplabel); rename("labels");
			labels = getImageID();
		}
		
		selectImage(id1); close();	
	}
	selectImage(id0); close();
	return labels;
}


function runHS(flt, args) {
	// Run a command on each frame of an hyperstack to correct a bug in ImageJ
	Stack.getDimensions(width, height, channels, slices, frames);
	src = getImageID();
	title = getTitle();
	for (n = 1; n <= frames; n++) {		
		selectImage(src);
		run("Duplicate...", "title=frame duplicate frames="+n);
		frame = getImageID();
		run(flt, args);
		if (n == 1) {
			rename("dst");
			dst = frame;			
		} else {
			run("Concatenate...", "title=dst open image1=dst image2=frame image3=[-- None --]");		
			dst = getImageID();			
		}
	}
	selectImage(src); close();
	selectImage(dst); rename(title);
	return dst;
}

function computeShell(labels) {
	/* compute ring masks if required
	 *  
	 */
	 
	selectImage(labels);	
	
	if (use_shell) {
		
		id1 = getImageID();
		
		run("Duplicate...", "title=min duplicate");
		//run("8-bit");
		runHS("Minimum 3D...", "x=1 y=1 z=1");
		id2 = getImageID();
		
		run("Duplicate...", "title=max duplicate");
		//run("8-bit");
		runHS("Maximum 3D...", "x=4 y=4 z=2");
		id3 = getImageID();
		
		imageCalculator("Subtract stack", id3, id2);
		selectImage(id1); close();
		selectImage(id2); close();
		selectImage(id3); rename("labels");
		setThreshold(0.5, 1.5);
		run("Convert to Mask", "background=Dark black");
		
	} else {
		
		runHS("Maximum 3D...", "x=2 y=2 z=1");	
		
	}	
	
	return getImageID();
}

function recordIntensity(labels, img) {
	// measurements in 3d of image in labels for each frame
	
	selectImage(img);
	Stack.getDimensions(width, height, channels, slices, frames);
	img_name = getTitle();
	
	selectImage(labels);
	labels_name = getTitle();
	
	Table.create("Result");
	
	for (frame = 1; frame <= frames; frame++) {
		
		selectImage(img);
		run("Duplicate...", "title=intensity duplicate channels=2 frames="+frame);
		
		selectImage(labels);
		run("Duplicate...", "title=mask duplicate frames="+frame);
		
		run("Intensity Measurements 2D/3D", "input=intensity labels=mask mean");
		selectWindow("intensity-intensity-measurements");
		value = Table.get("Mean", 0);		
		run("Close");
		
		selectWindow("intensity"); close();
		selectWindow("mask"); close();
		
		selectWindow("Result");
		Table.set("Frame", frame-1, frame);
		Table.set("Time [min]", frame-1, start_time + (frame + first_frame) * Stack.getFrameInterval()/60);
		Table.set("Mean intensity", frame-1, value);
		
	}
	Table.update;
}

function graphAndFit(name) {
	/*
	 * Create a graph and fit a model to the intensity recovery
	 *  - name: image name
	 *  - correct_bleaching: boolean which models
	 * Results
	 *  - plot
	 *  - table with parameters	 
	 */
	selectWindow("Result");
	x = Table.getColumn("Time [min]");
	y = Table.getColumn("Mean intensity");
		
	Plot.create("Intensity", "Time [min]", "Mean Intensity");
	Plot.setColor("#5555AA");	
	Plot.add("circle",x,y);
	Plot.setColor("#AA55AA");
	
	// initial guesses
	Array.getStatistics(y, ymin, ymax, mean, stdDev);
	a = (ymin + ymax) / 2;
	b = (ymax - ymin) / 2;
	c = (x[x.length-1] + x[0]) / 2;
	d = (x[x.length-1] - x[0]) / 10;
	e =  x[x.length-1] - x[0];
	
	if (correct_bleaching) {
		Fit.doFit("y=(a+b*Math.erf((x-c)/d))*exp(-x/e)", x, y,newArray(a,b,c,d,e));		
	} else {		
		Fit.doFit("Error Function", x, y, newArray(a,b,c,d));
	}
	a = Fit.p(0);
	b = Fit.p(1);
	c = Fit.p(2);
	d = Fit.p(3);
	R2 = Fit.rSquared;	
	
	fun = newArray(x.length);
	for (i = 0; i < fun.length; i++){
		fun[i] = Fit.f(x[i]);
	}
	Plot.add("line", x, fun);
	print("95% Recovery :" + 2.326 *d + " min");
	
	if (!isOpen("Summary")) {
		Table.create("Summary");
	} 
	selectWindow("Summary");
	nrow = Table.size;
	Table.set("Image", nrow, name);
	Table.set("Time mid point [min]", nrow, c);
	Table.set("Time constant [min]", nrow, d);
	Table.set("5%-95% recovery time [min]", nrow, 2.326 * d);
	Table.set("Intensity amplitude [au]", nrow, b);
	Table.set("Intensity offset [au]", nrow, a);
	if (correct_bleaching) {
		Table.set("Photo-bleaching [min]", nrow, e);	
	}
	Table.set("R squared", nrow, R2);
	Table.update;
}

function getDataSet(path) {
	if (matches(path, ".*image")) {	
		return getImageID();
	} else {
		run("Close All");
		open(path);		
		return getImageID();
	}
}

function saveAndFinih(path) {	
	if (!matches(path, ".*image")) {
		folder = File.getDirectory(path);
		fname = File.getNameWithoutExtension(path);		
		selectWindow("Result");
		Table.save(folder + File.separator + fname + ".csv");
		selectWindow("labels");
		saveAs("TIFF", folder + File.separator + fname + "-mask.tif");
	}
}

function main() {
	/* Entry point */
	setBatchMode("hide");	
	img = getDataSet(path);
	name = getTitle();
	labels = segmentAndTrackParasite(img);
	shell = computeShell(labels);
	recordIntensity(shell, img);
	graphAndFit(name);
	setBatchMode("exit and display");
	saveAndFinih(path);
	print("Done");
}

main();