//@Boolean(label="Use ring") use_ring
//@Float(label="start time [min]",value=23) start_time
//@Float(label="first frame [frame]",value=23) first_frame

/*
	Measure RNF213 coat fluorescence overtime
		
	Segment the parasites in channel 1 and measure mean
	intensity on a shell around the parasite in channel 2.
	
	We assume that there is only a parasite centered in the field of
	view at the start of the sequence.
	
	Installation:
	
	Open the macro using the script editor in Fiji.
	
	
	Usage:
	
	- Duplicate the region of interest with only 1 parasite
	- Run the macro
	- Save the table	
	
	
	
	Jerome Boulanger for Ana Crespillo Casado
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
	 * 	img : input image id
	 * 	
	 * Output:
	 * 	labels: labels window id
	 */
	Stack.getDimensions(width, height, channels, slices, frames);
	Stack.getUnits(X, Y, Z, Time, Value);
	getVoxelSize(dx, dy, dz, unit);
	print(width);
	run("Duplicate...", "title=binary duplicate channels=1");
	for (iter = 0; iter < 5; iter++) {
		run("Median 3D...", "x=2 y=2 z=1");
	}
	Stack.getStatistics(voxelCount, mean, min, max, stdDev);
	setThreshold(mean+stdDev, max+1);
	//setAutoThreshold("Default dark");	
	run("Convert to Mask", "background=Dark black");	
	id0 = getImageID();
	
	position = newArray(width / 2 * dx, height / 2 * dy, slices / 2 * dz);
	
	for (frame = 1; frame <= frames; frame++) {
		
		print("Processing frame " + frame);
		selectImage(id0);
		run("Duplicate...", "title=test duplicate frames="+frame);
		id1 = getImageID();
		
		run("Connected Components Labeling", "connectivity=6 type=[16 bits]");
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

function computeRing(labels) {
	// compute ring masks if required
	selectImage(labels);
	if (use_ring) {
		id1 = getImageID();
		run("Minimum 3D...", "x=1 y=1 z=1");
		run("Duplicate...", "duplicate");
		id2 = getImageID();
		run("Maximum 3D...", "x=4 y=4 z=2");
		imageCalculator("Subtract stack", id2, id1);
		selectImage(id1); close();
		rename("labels");
		return id2;	
	} else {
		run("Maximum 3D...", "x=2 y=2 z=1");
		return getImageID();
	}	
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

function graphAndFit() {
	selectWindow("Result");
	x = Table.getColumn("Time [min]");
	y = Table.getColumn("Mean intensity");
	// crop
	pos = Array.rankPositions(y);
	Array.print(pos);
	kmax = pos[pos.length-1];	
	Plot.create("Intensity", "Time [min]", "Mean Intensity");
	Plot.setColor("#5555AA");	
	Plot.add("circle",x,y);
	Plot.setColor("#AA55AA");
	Fit.doFit("Error Function", Array.trim(x,kmax), Array.trim(y,kmax));
	a = Fit.p(0);
	b = Fit.p(1);
	c = Fit.p(2);
	d = Fit.p(3);
	f = newArray(x.length);
	for (i=0;i<f.length;i++){
		f[i] = Fit.f(x[i]);
	}
	Plot.add("line", x, f);
	print("Time constant " + d + " min");
	if (!isOpen("Summary")) {
		Table.create("Summary");
	} 
	selectWindow("Summary");
	nrow = Table.size;
	Table.set("Time mid point [min]", nrow, c);
	Table.set("Time constant [min]", nrow, d);
	Table.set("Intensity amplitude [au]", nrow, b);
	Table.set("Intensity offset [au]", nrow, a);
	
	Table.update;
}


function main() {
	setBatchMode("hide");
	img = getImageID();
	labels = segmentAndTrackParasite(img);
	ring = computeRing(labels);
	recordIntensity(ring, img);
	graphAndFit();
	setBatchMode("exit and display");
	print("Done");
}

main();