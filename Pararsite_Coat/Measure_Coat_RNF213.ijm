//@File(label="Input",description="Input dataset tif or csv file, use 'image' to run on current image",value="image") path
//@Float(label="start time [min]",value=23) start_time
//@Float(label="first frame [frame]",value=23) first_frame
//@Boolean(label="Use shell") use_shell
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
	
	Processing an opened image:
	- Open the image to process
	- Duplicate the region of interest with only 1 parasite	
	- Run the macro with image as input
	- Save the table
	
	Processing a file
	- Open the image, crop the region to analyse and save them in a tif file
	- Keep track of the starting frame of the crop
	- Run the macro
	
	Processing several files
	- Open image, crop and export as tif, note in a csv file the first frame
	- The csv file should have the columns
	 Filename, Start time [min], First frame
	- Run the macro
	- Use the csv file as input
	
	
	
	
	Jerome Boulanger for Ana Crespillo Casado 2023
*/

function findNearestLabel(labels_frame, position) {
	/* Find nearest label and update position
	 * Return id of the image with the correct label
	 */
	selectImage(labels_frame);
	name = getTitle();
	Stack.getStatistics(voxelCount, mean, min, max, stdDev);
	if (max == 0) {print("Warning label not found"); return 0;}
	
	run("Analyze Regions 3D", "centroid surface_area_method=[Crofton (13 dirs.)] euler_connectivity=6");
	selectWindow(name+"-morpho");	
	xk = Table.getColumn("Centroid.X");
	yk = Table.getColumn("Centroid.Y");
	zk = Table.getColumn("Centroid.Z");
	labs = Table.getColumn("Label");	
	run("Close");
	// double check
	if (isOpen("test-labels-morpho")) {
		selectWindow("test-labels-morpho");
		run("Close");
	}
	
	dmin = 1e9;
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
	
	selectImage(labels_frame);
	run("Select Label(s)", "label(s)=" + labs[kstar]);
	
	selectWindow(name+"-keepLabels");
	labels_keep = getImageID();
	run("Remap Labels");
	
	//print("Found ", labs.length, "keep label", labs[kstar]);
	return labels_keep;
}


function runHS(flt, args, niter) {
	/* Run a command on each frame of an hyperstack to correct a bug in ImageJ */
	Stack.getDimensions(width, height, channels, slices, frames);
	src = getImageID();
	title = getTitle();
	for (n = 1; n <= frames; n++) {		
		selectImage(src);
		run("Duplicate...", "title=frame duplicate frames="+n);
		frame = getImageID();
		for (iter = 0; iter < niter; iter++) {
			run(flt, args);
		}
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

function watershed(mask){
	selectImage(mask);
	getVoxelSize(dx, dy, dz, unit);
	s = 0.5;
	tmask = getTitle();
	run("Chamfer Distance Map 3D", "distances=[Quasi-Euclidean (1,1.41,1.73)] output=[32 bits] normalize");
	selectWindow(tmask+"-dist");	
	run("Gaussian Blur 3D...", "x="+s/dx+" y="+s/dy+" z="+2*s/dz);
	run("Invert", "stack");
	run("Classic Watershed", "input="+tmask+"-dist mask="+tmask+" use min=0 max=255");
	selectWindow("watershed");
	rename(tmask+"-labels");
	res = getImageID();
	selectImage(tmask+"-dist"); close();
	return res;
}

function segmentAndTrackParasite(img) {
	/* Segment and track the parasite
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
	run("32-bit");
	run("Gaussian Blur 3D...", "x=0.7 y=0.7 z=0.7");
	for (i = 0; i < 5; i++) {				
		run("Median 3D...", "x=2 y=2 z=1");				
	}
	
	id0 = getImageID();
	
	Stack.getStatistics(voxelCount, mean, min, max, stdDev);
	setThreshold(mean + stdDev, max + 1);	
	run("Convert to Mask", "background=Dark black");	
	id0 = getImageID();
	
	position = newArray(width / 2 * dx, height / 2 * dy, slices / 2 * dz);
	
	for (frame = 1; frame <= frames; frame++) {
		
		selectImage(id0);
		run("Duplicate...", "title=test duplicate frames="+frame);
		id1 = getImageID();
		
		id2 = watershed(id1);	
		
		keeplabel = findNearestLabel(id2, position);
		selectImage(keeplabel);
		keepname = getTitle();
		
		selectImage(id2); close();
							
		if (frame > 1) {	
			run("Concatenate...", "title=labels open image1=labels image2="+keepname+" image3=[-- None --]");		
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

function computeShell(labels) {
	/* compute shell masks if required by use_shell input */
	t = 1.; // thickness
	s = 0.3; // shift
	
	a = 0.5 * (t - s);
	b = 0.5 * (t + s);	
		 
	selectImage(labels);	
	getVoxelSize(width, height, depth, unit);

	if (use_shell) {
		
		id1 = getImageID();
		
		run("Duplicate...", "title=min duplicate");		
		runHS("Minimum 3D...", "x="+Math.ceil(a/width)+" y="+Math.ceil(a/height)+" z="+Math.ceil(a/depth), 1);
		id2 = getImageID();
		
		run("Duplicate...", "title=max duplicate");		
		runHS("Maximum 3D...", "x="+Math.ceil(b/width)+" y="+Math.ceil(b/height)+" z="+Math.ceil(b/depth), 1);
		id3 = getImageID();
		
		imageCalculator("Subtract stack", id3, id2);
		selectImage(id1); close();
		selectImage(id2); close();
		selectImage(id3); rename("labels");
		setThreshold(0.5, 1.5);
		run("Convert to Mask", "background=Dark black");
		
	} else {
		
		runHS("Maximum 3D...", "x=2 y=2 z=1", 1);	
		
	}	
	setMinAndMax(0, 1);
	return getImageID();
}

function computeSumAtFrame(id, frame) {
	selectImage(id);
	Stack.getDimensions(width, height, channels, slices, frames);
	s = 0;
	run("Select None");
	for (slice = 1; slice <= slices; slice++) {		
		Stack.setPosition(1, slice, frame);
		s += getValue("RawIntDen");		
	}	
	return s;
}

function recordIntensity(labels, img, start_time, first_frame) {
	/* Measurements in 3d of image in labels for each frame */
	
	if (isOpen("test-labels-morpho")) {
		selectWindow("test-labels-morpho");		
		run("Close");
	}
	
	selectImage(img);
	Stack.getDimensions(width, height, channels, slices, frames);
	img_name = getTitle();
	
	selectImage(labels);
	labels_name = getTitle();
		
	selectImage(img);
	run("Duplicate...", "title=intensity duplicate channels=2");
	run("32-bit");
	img_prod = getImageID();	
	
	imageCalculator("Multiply stack", img_prod, labels);
	
	selectImage(img);
	run("Duplicate...", "title=intensity duplicate channels=2");
	run("32-bit");
	img_mask = getImageID();
	run("Median 3D...", "x=3 y=3 z=1");
	Stack.getStatistics(voxelCount, mean, min, max, stdDev);
	setThreshold(mean+stdDev, max+1);	
	run("Convert to Mask", "background=Dark black");	
	imageCalculator("Multiply stack", img_mask, labels);
			
	aframe = newArray(frames);
	atime = newArray(frames);
	amean = newArray(frames);
	afraction = newArray(frames);
	for (frame = 1; frame <= frames; frame++) {				
		s1 = computeSumAtFrame(img_prod, frame);
		s2 = computeSumAtFrame(labels, frame);
		s3 = computeSumAtFrame(img_mask, frame);		
		aframe[frame-1] = first_frame + frame - 1;
		atime[frame-1] = start_time + (frame - 1 + first_frame) * Stack.getFrameInterval()/60;
		if (s2 > 0) {
			amean[frame-1] = s1 / s2;
			afraction[frame-1] =  s3 / s2;		
		} else {
			amean[frame-1] = mean;
		}
	
	}
	Table.create("Intensities");
	selectWindow("Intensities");
	Table.setColumn("Frame", aframe);
	Table.setColumn("Time [min]", atime);
	Table.setColumn("Mean intensity", amean);
	Table.setColumn("Fraction", afraction);
	Table.update;
	selectImage(img_prod); close();
	selectImage(img_mask); close();
}

function recordIntensityOld(labels, img, start_time, first_frame) {
	/* Measurements in 3d of image in labels for each frame */
	
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
		Stack.getStatistics(voxelCount, mmean, mmin, mmax, mstdDev);
		
		if (mmax>0) {				
			run("Intensity Measurements 2D/3D", "input=intensity labels=mask mean");
			selectWindow("intensity-intensity-measurements");
			value = Table.get("Mean", 0);		
			run("Close");
		} else {
			print("Warning empty mask at frame " + frame);
			value = 0;
		}
		
		selectWindow("intensity"); close();
		selectWindow("mask"); close();
		
		selectWindow("Result");
		Table.set("Frame", frame-1, first_frame + frame - 1);
		Table.set("Time [min]", frame-1, start_time + (frame - 1 + first_frame) * Stack.getFrameInterval()/60);
		Table.set("Mean intensity", frame-1, value);		
	}
	Table.update;
}


function getObjectiveName() {
	txt = getMetadata("Info");
	lines = split(txt,"\n");
	for (i = 0; i < lines.length; i++) {
		if (matches(lines[i],"wsObjectiveName =.*")) {
			kv = split(lines[i],"=");		
			return kv[1];
		}		
	}
	return "unknown";
}

function getCameraName() {
	txt = getMetadata("Info");
	lines = split(txt,"\n");
	for (i = 0; i < lines.length; i++) {
		if (matches(lines[i],"wsCameraName =.*")) {
			kv = split(lines[i],"=");		
			return kv[1];
		}	
	}
	return "unknown";
}

function getNumericalAperture() {
	txt = getMetadata("Info");
	lines = split(txt,"\n");
	for (i = 0; i < lines.length; i++) {
		if (matches(lines[i],"Numerical Aperture =.*")) {
			kv = split(lines[i],"=");	
			return kv[1];
		}
	}
	return "unknown";
}


function graphAndFit(img, name, start_time, first_frame) {
	/* Create a graph and fit a model to the intensity recovery
	 *  - name: image name
	 *  - correct_bleaching: boolean which models
	 * Results
	 *  - plot
	 *  - table with parameters	 
	 */
	
	selectWindow("Intensities");
	x = Table.getColumn("Time [min]");
	y = Table.getColumn("Mean intensity");	
		
	Plot.create("Intensity", "Time [min]", "Mean Intensity");
	Plot.setColor("#5555AA");	
	Plot.add("circle",x,y);
	Plot.setColor("#AA55AA");
	
	// initial guesses
	Array.getStatistics(y, ymin, ymax, mean, stdDev);	
	a = (ymin + ymax) / 2;
	b = (ymax - ymin)/2;
	c = 0;	
	xstart = 0;
	xstop = 0;
	for (i = 0; i < x.length; i++) {
		if (y[i] > ymin + 0.1 * (ymax-ymin)  && xstart == 0) {
			xstart = x[i];
		}
		if (y[i] > a && c == 0) {
			c = x[i];
			print(i, x[i]);
		}
		if (y[i] > ymin + 0.9 * (ymax-ymin) && xstop == 0) {
			xstop = x[i];
		}
	}
	print("start",xstart, "end", xstop);
	if (xstop - xstart > 1 && xstop - xstart < (x[x.length-1] - x[0])/10) {
		d = (xstop - xstart) / 1.81;
	} else {
		d = (x[x.length-1] - x[0]) / 10;
	}
	print("c = ", c);
	print("d = ", d);
	e =  1 / (x[x.length-1] - x[0]);
	
	print("  initial values : [",a,b,c,d,e,"]");	
	// Error function model
	Fit.doFit("Error Function", x, y, newArray(a,b,c,d));
	a = Fit.p(0);
	b = Fit.p(1);
	c = Fit.p(2);
	d = Fit.p(3);	
	R2 = Fit.rSquared;	
	print("  erf fit: [",a,b,c,d,"] R2:", R2);
	// Error function + bleaching model
	if (correct_bleaching) {	
		Fit.doFit("y=a+b*Math.erf((x-c)/d)*exp(-e*x)", x, y, newArray(a,b,c,d,e));
		if (Fit.rSquared > R2) {
			a = Fit.p(0);
			b = Fit.p(1);
			c = Fit.p(2);
			d = Fit.p(3);
			e = Fit.p(4);
			R2 = Fit.rSquared;	
			print("  erf x bleach fit: [",a,b,c,d,e,"] R2:", R2);
		} else {
			e = 0;
			print("  erf x bleach fit failed R2:", Fit.rSquared);
		}
	}	
	
	fun = newArray(x.length);
	for (i = 0; i < fun.length; i++){
		fun[i] = Fit.f(x[i]);
	}
	Plot.add("line", x, fun);
	print("  > 95% Recovery :" + 2.326 *d + " min <");
	
	Table.setColumn("Model", fun);
	
	if (!isOpen("Summary")) {
		Table.create("Summary");
	} 
	selectWindow("Summary");
	nrow = Table.size;
	Table.set("Image", nrow, name);
	Table.set("Start time [min]", nrow,  start_time);
	Table.set("First frame", nrow,  first_frame);
	Table.set("Time mid point [min]", nrow, c);
	Table.set("Time constant [min]", nrow, d);
	Table.set("5%-95% recovery time [min]", nrow, 2.326 * d);
	Table.set("Intensity amplitude [au]", nrow, b);
	Table.set("Intensity offset [au]", nrow, a);
	if (correct_bleaching) {
		Table.set("Photo-bleaching [min]", nrow, 1/e);	
	}
	Table.set("R squared", nrow, R2);
	selectImage(img);
	Stack.getDimensions(width, height, channels, slices, frames);
	dt = Stack.getFrameInterval();
	getVoxelSize(dx, dy, dz, unit);
	camera = getCameraName();
	objective = getObjectiveName();
	na = getNumericalAperture();
	selectWindow("Summary");
	Table.set("Image width", nrow, width);
	Table.set("Image height", nrow, height);
	Table.set("Image slices", nrow, slices);
	Table.set("Image frames", nrow, frames);
	Table.set("Image pixel size", nrow, dx);	
	Table.set("Image spacing", nrow, dz);
	Table.set("Image frame interval", nrow, dt);
	Table.set("Objective", nrow, objective);
	Table.set("Camera", nrow, camera);
	Table.set("Numerical Aperture", nrow, na);
	Table.update;
}

function getDataSet(path) {
	/* Get the dataset given the path parameter */
	if (matches(path, ".*image")) {	
		print("Processing " + getTitle());
		return getImageID();
	} else {
		print("Processing " + path);
		run("Close All");
		open(path);		
		return getImageID();
	}
}

function saveAndFinish(img, shell, path) {	
	/* Save results if an image was opened */
	if (!matches(path, ".*image")) {
		folder = File.getDirectory(path);
		fname = File.getNameWithoutExtension(path);		
		selectWindow("Intensities");
		Table.save(folder + File.separator + fname + "-profile.csv");		
		selectImage(shell);
		run("16-bit");
		setMinAndMax(0, 1000);
		c3 = getTitle();
		selectImage(img);
		c = getTitle();
		run("16-bit");
		run("Split Channels");
		run("Merge Channels...","c1=[C1-"+c+"] c2=[C2-"+c+"] c3=["+c3+"]");
		saveAs("TIFF", folder + File.separator + fname + "-mask.tif");
	}
}

function processImage(path, start_time, first_frame) {
	/* Process the current image or a file */
	img = getDataSet(path);
	name = getTitle();
	print("- Tracking");
	labels = segmentAndTrackParasite(img);
	print("- Shell computation");
	shell = computeShell(labels);
	print("- Intensity measurement");
	recordIntensity(shell, img, start_time, first_frame);
	print("- Model fitting");
	graphAndFit(img, name, start_time, first_frame);
	saveAndFinish(img, shell, path);
}

function main() {
	/* Entry point */
	tic = getTime();
	setBatchMode("hide");	
	if (matches(path, ".*csv")) {
		print("Loading file list from", path);
		folder = File.getDirectory(path);
		Table.open(path);		
		filelist = Table.title();
		nrows = Table.size;
		for (row = 0; row < nrows; row++) {
			selectWindow(filelist);
			fname = Table.getString("Filename", row);
			start_time = Table.get("Start time [min]", row);
			first_frame = Table.get("First frame", row);
			print(fname, start_time, first_frame);
			processImage(folder + File.separator + fname, start_time, first_frame);
		}
	} else {
		processImage(path, start_time, first_frame);
	}
	setBatchMode("exit and display");
	toc = getTime();
	print("Done in " + (toc - tic) / 3600 + " seconds");
}

main();