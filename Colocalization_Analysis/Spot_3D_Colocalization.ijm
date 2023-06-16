// @File(label="Input",description="Use image to run on current image or test to run on a generated test image",value="image") path
// @Integer(label="Channel A", value=1) channela
// @Double(label="Spot size A [um]", value=0.2) spot_sizea
// @Double(label="Specificity A [log10]", value=1, style=slider, min=0, max=6) pfaa
// @Integer(label="Channel B", value=1) channelb
// @Double(label="Spot size B [um]", value=0.2) spot_sizeb
// @Double(label="Specificity B [log10]", value=1, style=slider, min=0, max=6) pfab
// @Double(label="Proximity threshold [um]",value=0.5) dmax
// @Boolean(label="Subpixel localization", value=true) subpixel
// @Boolean(label="Load coordinate from file", value=false) reload
// @Boolean(label="Close on exit", value=false) closeonexit
// @Boolean(label="Save coordinates", value=false) savecoordinates

/*
 * Detect 3D spots and count the number of pairs
 *
 *
 * Input:
 * - Channel A : first channel
 * - Channel B : second channel
 * - Spot size [um]: approximate size of the spot in micron
 * - Specificity [log10] : the higher the more strict (less detection)
 * - Proximity threshold [um] : the threshold to consider two points as close to each other
 * - subpixel localization : enable subpixel accuracy for localization
 * - load coordinates from file: if processing an image from disk, load the -points.csv file associated
 * - Close on exit: close the image and coordinates table when finished
 * - Save coordinates: save coordinates before exiting into a 'points.csv' file
 *
 *
 * If there is not opened image, a test image is generated
 *
 * The detection of the spots is done using a difference of Gaussian filters.
 * The localization of the spots is not subpixelic.
 * The voxel size is taken into account.
 *
 * Output:
 *
 * The results are saved in a table with the input parameters and
 * N(A) : number of spots in channel A
 * N(B) : number of spots in channel B
 * N(AB) : number of spots in B close to A
 * N(BA) : number of spots in A close to B
 * N(AB)/N(A) : ratio
 * N(BA)/N(B) : ratio
 *
 * Overlays are added to the stack to visualize the dections
 *
 * The macro doesn't provide any information whether or not the proximity count is
 * explained by randomness and the density of the points.
 *
 * Jerome for Anna and Andrei 2023
 *
 */

start_time = getTime();

mode = 0;
setBatchMode("hide");
if (matches(path,".*test")) {
	generateTestImage();
	channela = 1;
	channelb = 2;
	spot_size = 0.2;
	id = getImageID();
	basename = "test image";
	mode = 1;
} else if (matches(path,".*image")) {
	print("[Using opened image "+getTitle()+"]");
	basename = File.getNameWithoutExtension(getTitle);
	id = getImageID();
	mode = 2;
} else {
	basename = File.getNameWithoutExtension(path);
	if ((reload && !closeonexit) || !reload) {
		print("Loading image");
		print(path);
		//open(path);
		run("Bio-Formats Importer", "open=["+path+"] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
		id = getImageID();
	}
	mode = 3;
}

if (reload && mode!=3) {
	exit("Cannot reload coordinates, select an image file on disk.");
}

if (savecoordinates && mode!=3) {
	exit("Cannot save coordinates, select an image file on disk.");
}

if (closeonexit && mode != 3) {
	ok = getBoolean("Do you want to close the selected image "+getTitle+" on exit?");
	if (!ok) {
		closeonexit = false;
	}
}

tbl1 = "Spot3DColocalisation.csv";
tbl2 = basename + "-points.csv";

// define the name of the table
opath = File.getDirectory(path) + File.separator + File.getNameWithoutExtension(path) + "-points.csv";

// detect spots
if (reload && mode==3) {
	if (!File.exists(opath)) {
		exit("Could not find file" + opath);
	}
	A = loadCoordsTable(opath, channela);
	B = loadCoordsTable(opath, channelb);
} else {
	A = detect3DSpots(channela, spot_sizea, pfaa, subpixel);
	B = detect3DSpots(channelb, spot_sizeb, pfab, subpixel);
}

NA = A.length / 3;
NB = B.length / 3;

// find pairs of spots
NAB = colocalization(A, B, dmax);
NBA = colocalization(B, A, dmax);

// record results in table
recordInTable(tbl1, basename);

// records coordinates in table
coords2Table(tbl2, A, channela,spot_sizea);
coords2Table(tbl2, B, channelb,spot_sizeb);

// add overlays
if (!closeonexit) {
	Stack.getDimensions(width, height, channels, slices, frames);
	Overlay.remove();
	fontsize  = round(height / 20);
	setFont("Arial", fontsize);
	setColorFromLut(channela);
	coords2Overlay(A,spot_sizea,channela);
	Overlay.drawString("Channel "+channela+" N="+NA, 10, fontsize);
	setColorFromLut(channelb);
	coords2Overlay(B,spot_sizeb,channelb);
	Overlay.drawString("Channel "+channelb+" N="+NB, 10, 2*fontsize);
	Overlay.show();
	run("Select None");
	Stack.setSlice(round(slices/2));
	Stack.setChannel(channela);
	resetMinAndMax();
	Stack.setChannel(channelb);
	resetMinAndMax();
	active_channels = "";
	for (i=1;i<=channels;i++) {
		if (i==channela||i==channelb) {
			active_channels+="1";
		} else {
			active_channels+="0";
		}
	}
	Stack.setDisplayMode("composite");
	Stack.setActiveChannels(active_channels);
}

if (mode==3 && savecoordinates) {
	print("Saving file\n" + opath);
	selectWindow(tbl2);
	Table.save(opath);
}

if (closeonexit) {
	selectWindow(tbl2);
	run("Close");
	if (!reload) {
		selectImage(id);
		close();
	}
}

end_time = getTime();
run("Collect Garbage");
run("Collect Garbage");

print("Finished in " + (end_time-start_time)/1000 + " seconds.");
setBatchMode("exit and display");


function setColorFromLut(channel) {
	Stack.setChannel(channel);
	getLut(reds, greens, blues);
	setColor(reds[255],greens[255],blues[255]);
}

function recordInTable(tbl, name) {
	if (!isOpen(tbl)) {
		Table.create(tbl);
	}
	selectWindow(tbl);
	row = Table.size;
	Table.set("Name", row, name);
	Table.set("ChA", row, channela);
	Table.set("ChB", row, channelb);
	Table.set("Spot size A[um]", row, spot_sizea);
	Table.set("Specificity A", row, pfaa);
	Table.set("Spot size B [um]", row, spot_sizeb);
	Table.set("Specificity B", row, pfab);
	Table.set("Proximity threshold [um]", row, dmax);
	Table.set("N(A)", row, NA);
	Table.set("N(B)", row, NB);
	Table.set("N(AB)", row, NAB);
	Table.set("N(BA)", row, NBA);
	Table.set("N(AB)/NA", row, NAB/NA);
	Table.set("N(BA)/NB", row, NBA/NB);
	Table.update();
}

function coords2Table(tbl, coords, channel,size) {
	/*
	 * Records the coordinates in a table
	 *
	 * Parameter
	 *  tbl (string): table name
	 *  coords (array): (3,n) array with x,y,z in um
	 *  channel : spot channel
	 *  size : spot size used for detection
	 *
	 */
	if (!isOpen(tbl)) {
		Table.create(tbl);
	}
	selectWindow(tbl);
	N = coords.length / 3;	
	row = Table.size;
	for (i = 0; i < N; i++) {
		Table.set("Channel", row+i, channel);
		Table.set("X [um]", row+i, coords[3*i]);
		Table.set("Y [um]", row+i, coords[3*i+1]);
		Table.set("Z [um]", row+i, coords[3*i+2]);
		Table.set("Size [um]", row+i, size);
	}
	Table.update();
}

function loadCoordsTable(path, channel) {
	/*
	 * Load coordinates from a table
	 */
	Table.open(path);
	N = Table.size;
	coords = newArray(3*N);
	k = 0;
	for (i = 0; i < N; i++) {
		if (Table.get("Channel", i)==channel) {
			coords[3*k] = Table.get("X [um]", i);
			coords[3*k+1] = Table.get("Y [um]", i);
			coords[3*k+2] = Table.get("Z [um]", i);
			k++;
		}
	}
	return Array.trim(coords,3*k);
}

function colocalization(A,B,dmax) {
	/*
	 * Count number of proximity event
	 *
	 * A (array)   : coordinate of the 1st set of points (3,na)
	 * B (array)   : coordinate of the 2st set of points (3,nb)
	 * threshold (number) : distance threshold in physical units
	 *
	 * returns return number of close spots
	 */
	 ncoloc = 0;
	 for (i = 0; i < A.length/3; i++) {
		xi = A[3*i];
	 	yi = A[3*i+1];
	 	zi = A[3*i+2];
	 	dmin = 1e9;
	 	jmin = -1;
	 	for (j = 0; j < B.length/3; j++) {
	 		dx = (xi - B[3*j]);
	 		dy = (yi - B[3*j+1]);
	 		dz = (zi - B[3*j+2]);
	 		d = sqrt(dx * dx + dy * dy + dz * dz);
	 		if (d < dmin) {
	 			dmin = d;
	 			jmin = j;
	 		}

	 	}
	 	if (dmin < dmax) {
	 		ncoloc++;
	 	}
	 }
	 return ncoloc;
}

function coords2Overlay(coords,size,channel) {
	/*
	 * Add overlay o the image based on coordinates
	 *
	 * Parameter
	 * coords (array): (3,n) coordinates in [um]
	 * size : size of the circles in [um]
	 * channel : channel
	 */
	getVoxelSize(dx, dy, dz, unit);
	R = 2 * size / dx;
	for (i = 0; i < coords.length / 3; i++) {
		x = coords[3*i] / dx;
		y = coords[3*i+1] / dy;
		z = coords[3*i+2] / dz;
		Stack.setSlice(round(z));
		Overlay.drawEllipse(x-R, y-R, 2*R+1, 2*R+1);
		Overlay.setPosition(0, z, 0);
		Overlay.add;
	}
	Overlay.show();
}

function detect3DSpots(channel, size, pfa, subpixel) {
	/*
	 * detect spots and return coordinates in physical units as an array
	 * 
	 * Parameter
	 *  channel (number) : index of the channel
	 *  pfa (number): probability of false alarm in -log10
	 *  subpixel (boolean): use subpixel localization
	 *  
	 *  Returns
	 *   an array (3,n) of coordinates (x,y,z) in physical units
	 */

	run("Select None");
	getVoxelSize(dx, dy, dz, unit);

	sigma1 = size;
	sigma2 = 2 * size;	
	id0 = getImageID;

	// compute a difference of Gaussian
	run("Duplicate...", "title=id1 duplicate channels="+channel);
	run("32-bit");
	id1 = getImageID;
	//run("Square Root","stack");
	run("Gaussian Blur 3D...", "x="+sigma1/dx+" y="+sigma1/dy+" z="+2.5*sigma1/dz);

	run("Duplicate...", "title=id2 duplicate");
	id2 = getImageID;
	run("Gaussian Blur 3D...", "x="+sigma2/dx+" y="+sigma2/dy+" z="+2.5*sigma2/dz);
	imageCalculator("Subtract 32-bit stack",id1,id2);
	selectImage(id2);close();

	// local maxima
	selectImage(id1);
	run("Duplicate...", "title=id3 duplicate");
	id3 = getImageID;
	run("Maximum 3D...", "x=1 y=1 z=1");
	imageCalculator("Subtract 32-bit stack",id3,id1);
	run("Macro...", "code=v=(v==0) stack");

	// threshold
	selectImage(id1);
	if (nSlices>1) {
		Stack.getStatistics(voxelCount, mean, min, max, stdDev);
	} else {
		getStatistics(area, mean, min, max, stdDev, histogram);
	}
	lambda = -2*normppf(pow(10,-pfa));
	threshold = mean + lambda * stdDev;	
	//print("specicity: "+pfa+" lambda: "+lambda+" threshold: "+threshold);
	run("Macro...", "code=v=(v>="+threshold+") stack");

	// combine local max and threshold
	imageCalculator("Multiply stack",id1,id3);
	selectImage(id3);close();
	setThreshold(0.5, 1);
	run("Make Binary", "method=Default background=Dark black");

	// count
	if (nSlices>1) {
		Stack.getStatistics(count, mean, min, max, stdDev);
	} else {
		getStatistics(count, mean, min, max, stdDev, histogram);
	}
	N =  mean * count / max;

	// localize
	Stack.getDimensions(width, height, channels, slices, frames);
	coords = newArray(3*N);
	j = 0;
	for (slice = 1; slice <= slices; slice++) {
		Stack.setSlice(slice);
		setThreshold(128,256);
		run("Create Selection");
		if (Roi.size != 0) {
			Roi.getContainedPoints(xpoints, ypoints);
			for (i = 0; i < xpoints.length; i++) {
				coords[3*j] = xpoints[i];
				coords[3*j+1] = ypoints[i];
				coords[3*j+2] = slice;
				j = j + 1;
			}
		}
	}
	coords = Array.trim(coords,3*j);
	close();

	// subpixel localization and conversion to physical units
	selectImage(id0);
	Stack.setChannel(channel);
	delta = newArray(dx,dy,dz);
	N = coords.length / 3;
	for (j = 0; j < N; j++) {
		if (subpixel) {
			p = blockCoordinateGaussianFit(Array.slice(coords,3*j,3*j+3));
		} else {
			p = Array.slice(coords,3*j,3*j+3);
		}
		for (k = 0; k < 3; k++) {
			coords[3 * j + k] = p[k] * delta[k];
		}
	}

	return coords;
}

function normcdf(x) {
	/*
	 * Normal umulative density function approximation
	 */
	s = 0;
	n = 1;
	for (i = 0; i < 100; i++) {
		n *= (2*i+1);
		s += pow(x,2*i+1) / n;
	}
	return 0.5 + 1. / sqrt(2.*PI) * exp(-0.5*x*x) * (s);
}

function normppf(p) {
	/*
	 * Percent point function (norm inv)
	 * Newton's method
	 *
	 * Parameter
	 * p : probability
	 *
	 * Returns
	 * the value x such that cdf(x) < p
	 */
	x = 0;
	for (i = 0; i < 50; i++) {
		delta = (normcdf(x) - p) / exp(-0.5*x*x) * sqrt(2 * PI);
		x = x - delta;
		if (abs(delta) < 1e-6){
			break;
		}
	}
	return x;
}


function blockCoordinateGaussianFit(p) {
	/*
	 * Fit a the location of a Gaussian spot
	 *
	 * Parameters:
	 * p (array): (3,) array with (x,y,z) coordinates
	 *
	 * Returns:
	 * p (array): (3,) array with (x,y,z) coordinates
	 *
	 */
	Stack.getDimensions(width, height, channels, slices, frames);
	n = 5;
	x = p[0];
	y = p[1];
	z = p[2];
	for (iter = 0; iter < 5; iter++) {
		a = getLineIntensity(x,y,z,n,0,0,n);
		Fit.doFit("Gaussian", Array.slice(a,0,n), Array.slice(a,3*n,4*n));
		x = Fit.p(2);
		a = getLineIntensity(x,y,z,0,n,0,n);
		Fit.doFit("Gaussian", Array.slice(a,n,2*n), Array.slice(a,3*n,4*n));
		y = Fit.p(2);
		if (slices > n) {
			a = getLineIntensity(x,y,z,0,0,n,n);
			Fit.doFit("Gaussian", Array.slice(a,2*n,3*n), Array.slice(a,3*n,4*n));
			z = Fit.p(2);
		}
	}
	if (abs(x-p[0]) < 1 && abs(y-p[1]) < 1 && abs(z-p[2]) < 1) {
		return newArray(x,y,z);
	} else {
		return p;
	}
}


function getLineIntensity(x,y,z,u,v,w,n) {
	/*
	 * Measure intensity along a line
	 *
	 * Parameters
	 * x : x-coordinates
	 * y : y-coordinates
	 * z : z-coordinates
	 * u : x-value of line vector
	 * v : y-value of line vector
	 * w : z-value of line vector
	 * n : number of sampling points
	 *
	 * Returns
	 * array (n,4) with x,y,z,intensity
	 * x  = Array.slice(0,n) etc..
	 */
	a = newArray(n);
	for (i = 0; i < n; i++) {
		xi = x + (i/(n-1)-0.5) * u;
		yi = y + (i/(n-1)-0.5) * v;
		zi = z + (i/(n-1)-0.5) * w;
		xi = minOf(width-1, maxOf(xi,0));
		yi = minOf(height-1, maxOf(yi,0));
		zi = minOf(slices, maxOf(zi,1));
		Stack.setSlice(zi);
		intensity = getPixel(xi, yi);
		a[i] = xi;
		a[n+i] = yi;
		a[2*n+i] = zi;
		a[3*n+i] = intensity;
	}
	return a;
}

function generateTestImage() {
	w = 200;
	h = 200;
	d = 100;
	n = 100;
	newImage("Test image", "16-bit grayscale-mode", w, h, 2, d, 1);
	run("Macro...", "code=100 stack");
	N1 = 0;
	N2 = 0;
	N12 = 0;
	b = 10;
	for (i = 0; i < n; i++) {
		x = b + round((w-2*b) * random);
		y = b + round((h-2*b) * random);
		z = b + round((d-2*b) * random);
		a_on = false;
		b_on = false;
		if (random > 0.2) {
			Stack.setPosition(1, z, 1);
			setPixel(x,y,10000);
			N1++;
			a_on = true;
		}
		if (random > 0.2) {
			Stack.setPosition(2, z, 1);
			setPixel(x,y,10000);
			N2++;
			b_on = true;
		}
		if (a_on && b_on) {
			N12++;
		}
	}
	run("Maximum 3D...", "x=1 y=1 z=1");
	sigma = 0.2;
	dxy = 0.11;
	dz  = 0.33;
	run("Gaussian Blur 3D...", "x="+sigma/dxy+" y="+sigma/dxy+" z="+3*sigma/dz);
	run("Add Noise", "stack");
	run("Properties...", "channels=2 slices=100 frames=1 pixel_width="+dxy+" pixel_height="+dxy+" voxel_depth="+dz);
	print("[Simulated]\nNA : " + N1 + "  NB : " + N2 + " NAB : " + N12);
	Stack.setChannel(1);
	setMinAndMax(50, 500);
	Stack.setChannel(2);
	setMinAndMax(50, 500);
}
