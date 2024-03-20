// @File(label="Input",description="Use image to run on current image or test to run on a generated test image",value="image") path
// @String(label="Channels", value="1,2,3",description="comma separated list of channels indices") channels_str
// @String(label="Spot Size", value="0.5,0.5,0.5",description="comma separated list of spot size in microns") spot_size_str
// @String(label="Feature", choices={"DoG","LoG","Top hat"}) feature
// @String(label="Specificity[-log10]", value="3,3,3",description="comma separated list of specificity (>0)") specificity_str
// @String(label="Mask channel",value="none") mask_str
// @Double(label="Proximity threshold [um]",value=0.5,description="define the max distance between two nearby spots") dmax
// @Boolean(label="Subpixel localization", value=true, description="enable subpixel localization") subpixel
// @Boolean(label="Close on exit", value=false,description="close after processing") closeonexit
// @Boolean(label="Save coordinates", value=false) savecoordinates
// @Boolean(label="Z project", value=false) dozproject
/*
 * Multiple channel spot 3d detection and colocalization
 *
 * Coordinate based colocalization based on a distance threshold.
 *
 * The input can be an image stack or a csv file with (c,x,y,z) formats.
 * use the word 'image' for processing the current image or 'test' to run a test as input.
 *
 * Channels are comma separated channel indices
 *
 * A channel with intensity > 0 can be used as a mask.
 *
 * The reported counts are the number of localization corresponding to each combination of channels.
 *
 * To get colocalization ratios, one need to compute for example Ch1Ch2 / (Ch1Ch2 + Ch1)
 *
 */

function getMode() {
	if (matches(path,".*test")) {
		path = "Test Image.tif";
		channels_str = "1,2,3";
		spot_size_str = "0.2,0.2,0.2,0.2";
		specificity_str = "1,1,1,1";
		mask_str = "";
		return 1;
	} else if (matches(path,".*image")) {
		if (nImages==0) {
			Dialog.create("Please select a file");
			Dialog.addFile("Input file", File.getDefaultDir);
			Dialog.show();
			path = Dialog.getString();
			open(path);
		}
		path = getTitle();
		return 2;
	} else if (matches(path,".*csv")) {
		return 3;
	} else {
		return 4;
	}
}

function parseCSVInt(csv) {
	a = split(csv,",");
	for (i = 0; i < a.length; i++) {
		a[i] = parseInt(a[i]);
	}
	return a;
}

function parseCSVFloat(csv) {
	a = split(csv,",");
	for (i = 0; i < a.length; i++) {
		a[i] = parseFloat(a[i]);
	}
	return a;
}

function parseMaskStr(str) {
	/* parse mask input str and returns nan is there is no channel*/
	return parseInt(str);
}

function getMask(channel) {
	if (isNaN(channel)) { return 0;	}
	print("Using channel " + channel + " as segmentation mask.");
	run("Duplicate...", "title=segmentation duplicate channels="+channel);
	run("Macro...", "code=v=(v>=1) stack");
	run("32-bit");
	return getImageID();
}

function generateTestImage(channels) {
	/* generate a test image with 3d spots and c channels */
	run("Close All");
	w = 200;
	h = 200;
	d = 100;
	n = 100;
	print("Generate " + n + " candidate locations.");
	c = channels.length;
	name = File.getNameWithoutExtension(path);
	newImage(name, "16-bit grayscale-mode", w, h, c, d, 1);
	run("Macro...", "code=100 stack");
	b = 10; // border
	codes = newArray(n*c*c); // 1 code per detection
	correct = newArray(n*c); // 1 code per location
	sets = powerSet(c);
	agg0 = newArray(sets.length);
	k = 0; // detection counter
	for (i = 0; i < n; i++) {
		x = b + round((w-2*b) * random);
		y = b + round((h-2*b) * random);
		z = b + round((d-2*b) * random);
		values = newArray(c);

		for (j = 0; j < c; j++) {
			if (random > 0.2) {
				Stack.setPosition(j+1, z, 1);
				setPixel(x,y,10000);
				values[j] = 1;
			}
		}

		// write an entry for each positive value with all the codes
		for (j0 = 0; j0 < c; j0++) if (values[j0]==1) {
			for (j = 0; j < c; j++) {
				codes[j+c*k] = values[j];
			}
			k = k + 1;
		}

		// write an entry for each location
		for (j = 0; j < c; j++) {
			correct[j+c*i] = values[j];
		}

		// compute directly the counts for each code
		indx = parseInt(array2str(values),2);
		agg0[indx] = agg0[indx] + 1;

	}
	codes = Array.trim(codes, c*k);
	run("Maximum 3D...", "x=1 y=1 z=1");
	sigma = 0.2;
	dxy = 0.11;
	dz  = 0.33;
	run("Gaussian Blur 3D...", "x="+sigma/dxy+" y="+sigma/dxy+" z="+3*sigma/dz);
	run("Add Noise", "stack");
	run("Properties...", "channels="+c+" slices="+d+" frames=1 pixel_width="+dxy+" pixel_height="+dxy+" voxel_depth="+dz);
	//print("[Simulated]\nNA : " + N1 + "  NB : " + N2 + " NAB : " + N12);
	Stack.setSlice(round(d/2));
	for (i = 1; i <= c; i++) {
		Stack.setChannel(i);
		setMinAndMax(50, 1000);
	}
	// showCodes("gt codes", codes, channels);
	sets = powerSet(c);
	counts = countsCodeBySet(codes, sets, channels);
	// showCounts("gt counts", counts, sets, channels);
	agg = aggregateCountsPerSet(counts, sets, channels);
	/*
	appendRecord(agg, sets, channels, "ground truth", newArray(c), newArray(c));
	Array.print(agg);

	counts2 = countsCodeBySet(correct, sets, channels);
	agg2 = aggregateCountsPerSet(counts2, sets, channels);
	appendRecord(agg2, sets, channels, "ground truth*", newArray(c), newArray(c));
	*/

	// compute ground truth from localization
	agg3 = correctAggCounts(agg, sets, channels);
	agg4 = divideBySetCardinality(agg3,sets,channels);
	appendRecord(agg4, sets, channels, "ground truth", newArray(c), newArray(c));

	// actual ground truth
	appendRecord(agg0, sets, channels, "ground truth*", newArray(c), newArray(c));
	//agg5 = correctAggCounts(agg2, sets, channels);
	//appendRecord(agg4, sets, channels, "ground truth***", newArray(c), newArray(c));
}

function array2str(a) {
	/* Convert an array to a string */
	s = "";
	for (i = 0; i < a.length; i++) {
		s = s + toString(a[i]);
	}
	return s;
}

function printAgg(agg, sets, channels) {
	/* Print the counts per set*/
	c = channels.length;
	for (k = 0; k < pow(2, c); k++) {
		code = Array.slice(sets,k*c,(k+1)*c);
		print("" + k + " " + array2str(code) + ":" + agg[k]);
	}
}



function testIntersection(code1, code2) {
	/* Test if two code intersect */
	if (code2.length==0) return false;
	for (i = 0; i < code1.length; i++) {
		if (code2[i] == 1 && code1[i] != 1) {
			return false;
		}
	}
	return true;
}

function correctAggCounts(agg, sets, channels) {
	/* Correct aggregated results */
	c = channels.length;
	agg2 = Array.copy(agg);
	// for each code i:a
	for (j = pow(2, c)-1; j > 1; j--) {
		a = Array.slice(sets, j*c, (j+1)*c);
		//print("\n " + j + " : " + array2str(a) + "");
		// for each other code j:b
		for (i = j-1; i > 0; i--) {
			b = Array.slice(sets, i*c, (i+1)*c);
			t = testIntersection(a, b);
			//print("    " + i + " : " + array2str(b) + " ~ " + t);
			if (t) {
				agg2[i] = agg2[i] - agg2[j];
			}
		}
	}
	return agg2;
}

function divideBySetCardinality(agg,sets,channels) {
	/* Divide aggregated counts by cardinality of the set */
	agg2 = Array.copy(agg);
	c = channels.length;
	for (k = 0; k < agg.length; k++) {
		code = Array.slice(sets, k*c, (k+1)*c);
		card = 0;
		for (i = 0; i < code.length; i++) {
			card += code[i];
		}
		agg2[k] = agg[k] / card;
	}
	return agg2;
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
	if (Math.log10(p) < -10) {
		return -sqrt(-4.2*Math.log10(p)+4.5);
	}
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

function normcdf(x) {
	/* Normal cumulative density function approximation	 */
	s = 0;
	n = 1;
	for (i = 0; i < 100; i++) {
		n *= (2*i+1);
		s += pow(x,2*i+1) / n;
	}
	return 0.5 + 1. / sqrt(2.*PI) * exp(-0.5*x*x) * (s);
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
	/* Measure intensity along a line
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

function addBorder(id, border, value) {
	/* Add borders to the image
	 *
	 */
	Stack.getDimensions(width, height, channels, slices, frames);
	run("Canvas Size...", "width="+(width+2*border)+" height="+(height+2*border)+" position=Center");
	run("Macro...", "code=v=(v==0)*100+(v!=0)*v stack");

	for (k = 0; k < border; k++) {
		Stack.setSlice(1);
		run("Add Slice", "add=slice prepend");
		Stack.setSlice(nSlices);
		run("Add Slice", "add=slice");
	}

	setColor(value);
	for (k = 1; k <= border; k++) {
		for (c = 1; c <= channels; c++) {
			Stack.setSlice(k);
			Stack.setChannel(c);
			fill();
			Stack.setSlice(border+slices+k);
			Stack.setChannel(c);
			fill();
		}
	}

	id0 = getImageID();

	run("Duplicate...", "duplicate");
	run("Gaussian Blur 3D...", "x="+2*border+" y="+2*border+" z="+2*border);
	id1 = getImageID();

	// copy the inner part into the new image
	for (c = 1; c <= channels; c++) {
		for (k = 1; k <= slices; k++) {

			selectImage(id0);
			Stack.setSlice(k+border);
			Stack.setChannel(c);
			makeRectangle(border, border, width, height);
			run("Copy");

			selectImage(id1);
			Stack.setSlice(k+border);
			Stack.setChannel(c);
			makeRectangle(border, border, width, height);
			run("Paste");
		}
	}
	selectImage(id0); close();
	selectImage(id1);
	run("Select None");
	return id1;
}

function cropBorder(id, border) {
	Stack.getDimensions(width, height, channels, slices, frames);
	makeRectangle(border, border, width-2*border, height-2*border);
	run("Duplicate...", "duplicate slices="+(border+1)+"-"+(slices-border));
	run("Select None");
	return getImageID();
}

function blobDOG(channel, square_root, size) {
	/* Compute a difference of Gaussian and returns the id */
	
	getVoxelSize(dx, dy, dz, unit);
	sigma1 = size;
	sigma2 = 3 * size;
	print("    scale: " + size + "um, [" + sigma1/dx + ","+ sigma1/dx+ ","+ 2.5*sigma1/dz+"] px");
	run("Duplicate...", "title=id1 duplicate channels="+channel);
	run("32-bit");
	id1 = getImageID();
	if (square_root) {
		run("Square Root","stack");
	}
	run("Gaussian Blur 3D...", "x="+sigma1/dx+" y="+sigma1/dy+" z="+2.5*sigma1/dz);

	run("Duplicate...", "title=id2 duplicate");
	id2 = getImageID;
	run("Gaussian Blur 3D...", "x="+sigma2/dx+" y="+sigma2/dy+" z="+2.5*sigma2/dz);
	imageCalculator("Subtract 32-bit stack",id1,id2);
	selectImage(id2);close();
	selectImage(id1);
	run("Top Hat...", "radius=10 stack");

	// naive border correction
	Stack.setSlice(1);
	run("Multiply...", "value=0.5");
	Stack.setSlice(nSlices);
	run("Multiply...", "value=0.5");
	return id1;
}

function blobLOG(channel, square_root, size) {
	/* Compute a laplacian of Gaussian and returns the id */
	sigma1 = size;
	getVoxelSize(dx, dy, dz, unit);
	print("    scale: " + size + "um, [" + sigma1/dx + ","+ sigma1/dx+ ","+ 2.5*sigma1/dz+"] px");
	run("Duplicate...", "title=id1 duplicate channels="+channel);
	run("32-bit");
	id1 = getImageID();
	if (square_root) {
		run("Square Root","stack");
	}
	run("Gaussian Blur 3D...", "x="+sigma1/dx+" y="+sigma1/dy+" z="+2.5*sigma1/dz);
	run("FeatureJ Laplacian", "compute smoothing=0.75");
	selectWindow("id1 Laplacian");
	id2 = getImageID();
	selectImage(id1); close();
	selectImage(id2);
	run("Multiply...", "value=-1 stack");
	rename("id1");

	// naive border correction
	Stack.setSlice(1);
	run("Multiply...", "value=0.5");
	Stack.setSlice(nSlices);
	run("Multiply...", "value=0.5");
	return id2;
}

function blobTH(channel, square_root, size) {
	/* Top hat detector */
	getVoxelSize(dx, dy, dz, unit);
	sx = maxOf(3 * size / dx, 2);
	sy = maxOf(3 * size / dy, 2);
	sz = maxOf(3 * size / dz, 2);
	run("Duplicate...", "title=id1 duplicate channels="+channel);
	run("32-bit");
	id1 = getImageID();
	if (square_root) {
		run("Square Root","stack");
	}
	run("Gaussian Blur 3D...", "x=0.75 y=0.75 z=0.75");
	run("Duplicate...", "duplicate");
	run("Minimum 3D...", "x="+sx+" y="+sy+" z="+sz);
	id2 = getImageID();
	imageCalculator("Subtract stack", id1, id2);
	selectImage(id2); close();
	selectImage(id1);
	return id1;
}

function localMaxima3D(id, size) {
	/* Compute local maxima and return id */
	selectImage(id);
	getVoxelSize(dx, dy, dz, unit);
	sx = maxOf(2 * size / dx, 1);
	sy = maxOf(2 * size / dy, 1);
	sz = maxOf(2 * size / dz, 1);
	run("Duplicate...", "title=id3 duplicate");
	dst = getImageID;
	run("Maximum 3D...", "x="+sx+" y="+sy+" z="+sz);
	imageCalculator("Subtract 32-bit stack", dst, id);
	run("Macro...", "code=v=(v==0) stack");
	//Stack.setSlice(1); setColor(0); fill();
	//Stack.setSlice(nSlices); setColor(0); fill();
	return dst;
}

function thresholdAuto(id, pfa) {
	/* threshold image id in place */
	selectImage(id);
	run("Select None");
	if (nSlices > 1) {
		Stack.getStatistics(voxelCount, mean, min, max, stdDev);
	} else {
		getStatistics(area, mean, min, max, stdDev, histogram);
	}

	lambda = -2 * normppf(pow(10, -pfa));
	threshold = mean + lambda * stdDev;

	print("    moments: [" + mean + ", " + stdDev + "]");
	print("    specificity: " + pfa + ", quantile: "+lambda+", threshold: "+threshold);
	run("Macro...", "code=v=(v>="+threshold+") stack");
}

function localizeSpots(id) {
	/* localize spots and return coordinates as a (n,4) array */

	setThreshold(0.5, 1);
	run("Make Binary", "method=Default background=Dark black");

	if (nSlices>1) {
		Stack.getStatistics(count, mean, min, max, stdDev);
	} else {
		getStatistics(count, mean, min, max, stdDev, histogram);
	}

	N =  mean * count / max;

	Stack.getDimensions(width, height, channels, slices, frames);
	coords = newArray(4*N);
	j = 0;
	for (slice = 1; slice <= slices; slice++) {
		Stack.setSlice(slice);
		setThreshold(128,256);
		run("Create Selection");
		if (Roi.size != 0) {
			Roi.getContainedPoints(xpoints, ypoints);
			for (i = 0; i < xpoints.length; i++) {
				coords[4*j]   = channel_idx;
				coords[4*j+1] = xpoints[i];
				coords[4*j+2] = ypoints[i];
				coords[4*j+3] = slice;
				j = j + 1;
			}
		}
	}
	coords = Array.trim(coords,4*j);
	return coords;
}

function refineLocalization(id, channel, coords) {
	selectImage(id);
	getVoxelSize(dx, dy, dz, unit);
	Stack.setChannel(channel);
	delta = newArray(dx,dy,dz);
	N = coords.length / 4;
	for (j = 0; j < N; j++) {
		if (subpixel) {
			p = blockCoordinateGaussianFit(Array.slice(coords,4*j+1,4*j+4));
		} else {
			p = Array.slice(coords,4*j+1,4*j+4);
		}
		for (k = 0; k < 3; k++) {
			coords[4 * j + k + 1] = p[k] * delta[k];
		}
	}
}

function detect3DSpots(channels_list, feature, size_list, pfa_list, channel_idx, subpixel, mask_id) {
	/*
	 * detect spots and return coordinates in physical units as an array
	 *
	 * Parameter
	 *  channels_list (array)  : list of indices of the channel in the image stack
	 *  size_list (array)
	 *  pfa_list (array) : probability of false alarm in -log10
	 *  channel_idx (int) : index of the channel in the array
	 *  subpixel (boolean): use subpixel localization
	 *  mask_channel (int or NaN): index of the channel
	 *
	 *  Returns
	 *   an array (4,n) of coordinates (c,x,y,z) in physical units
	 */

	channel = channels_list[channel_idx];
	size = size_list[channel_idx];
	pfa = pfa_list[channel_idx];

	print(" - detect spots in channel " + channel);

	run("Select None");

	id0 = getImageID;

	// compute a difference of Gaussian
	if (feature == "DOG") {
		id1 = blobDOG(channel, true, size);
	} else if (feature == "LOG") {
		id1 = blobLOG(channel, true, size);
	} else if (feature == "Top hat") {
		id1 = blobTH(channel, true, size);
	} else {
		run("Duplicate...", "title=id1 duplicate channels="+channel);
		id1 = getImageID();
	}

	// local maxima
	id3 = localMaxima3D(id1, size);

	// combine result with mask
	if (mask_id < 0) {
		imageCalculator("Multiply stack", id3, mask_id);
	}

	// threshold
	thresholdAuto(id1, pfa);

	// combine local max and threshold
	imageCalculator("Multiply stack", id1, id3);
	selectImage(id3); close();

	coords = localizeSpots(id1);

	refineLocalization(id1, channel, coords);

	selectImage(id1); close();

	print("   number of spots: " + coords.length/4);
	return coords;
}

function coords2Table(tbl, coords, channels, sizes) {
	/*
	 * Records the coordinates in a table
	 *
	 * Parameter
	 *  tbl (string): table name
	 *  coords (array): (3,n) array with x,y,z in um
	 *  channels (array): spot channel list
	 *  size (number): spot size used for detection
	 *
	 */

	Table.create(tbl);
	selectWindow(tbl);
	N = coords.length / 4;
	row = Table.size;
	for (i = 0; i < N; i++) {
		Table.set("Channel",   row+i, channels[coords[4*i]]);
		Table.set("X [um]",    row+i, coords[4*i+1]);
		Table.set("Y [um]",    row+i, coords[4*i+2]);
		Table.set("Z [um]",    row+i, coords[4*i+3]);
		Table.set("Size [um]", row+i, sizes[coords[4*i]]);
	}
	Table.update();
}

function loadCoordsTable(path, channels) {
	/*
	 * Load coordinates from a table
	 *
	 * returns
	 * a (4,n) array with the channel index in the channels
	 */
	Table.open(path);
	N = Table.size;
	coords = newArray(4*N);
	for (i = 0; i < N; i++) {
		channel = Table.get("Channel", i);
		for (idx = 0; idx < channels.length; idx++) {
			if (channels[idx] == channel) { break; }
		}
		coords[4*i]   = idx;
		coords[4*i+1] = Table.get("X [um]", i);
		coords[4*i+2] = Table.get("Y [um]", i);
		coords[4*i+3] = Table.get("Z [um]", i);
	}
	return coords;
}

function detectSpotsInAllChannels(id, channels, feature, specificity, spot_size, mask_channel) {
	setBatchMode("hide");
	selectImage(id);
	print("Detect 3D spots in image " + getTitle);
	coords = newArray(0);
	for (idx = 0; idx < channels.length; idx++) {
		current_coords = detect3DSpots(channels, feature, spot_size, specificity, idx, subpixel, mask_channel);
		coords = Array.concat(coords, current_coords);
	}
	setBatchMode("exit and display");
	return coords;
}


function loadData(mode, channels, feature, specificity, spot_size, mask_channel) {
	/* Load points data and set the basename */
	if (mode==1) {
		print("Test image");
		generateTestImage(channels);
		id = getImageID();
		basename = "test image";
		coords = detectSpotsInAllChannels(id, channels, feature, specificity, spot_size, mask_channel);
		//if (isOpen(basename+"-points.csv")) {selectWindow(basename+"-points.csv");run("Close");}
		//coords2Table(basename+"-points.csv", coords, channels, spot_size);
	} else if (mode==2) {
		print("[Using active image "+getTitle()+"]");
		basename = File.getNameWithoutExtension(getTitle);
		id = getImageID();
		mask_id = getMask(mask_channel);
		coords = detectSpotsInAllChannels(id, channels, feature, specificity, spot_size, mask_id);
		if (mask_id < 0) {
			selectImage(mask_id); close();
		}
		//coords2Table(basename+"points.csv", coords, channels, spot_size);
	} else if (mode==3) {
		basename = File.getNameWithoutExtension(path);
		print("Loading csv " + path + " table");
		open(path);
		coords = loadCoordsTable(path);
	} else {
		basename = File.getNameWithoutExtension(path);
		print("Opening image");
		print(path);
		run("Bio-Formats Importer", "open=["+path+"] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
		id = getImageID();
		coords = detectSpotsInAllChannels(id, channels, feature, specificity, spot_size);
		//coords2Table(basename+"-points.csv", coords, channels, spot_size);
	}
	return coords;
}


function computeCodes(coords, channels, dmax) {
	/* Count number of close by spots for each spots in every channel
	 *
	 * Parameters
	 *  coords (array)   : (4,N) coordinates c,x,y,z
	 *  channels (array) : indices of the channels
	 *  dmax (float)     : distance threshold
	 * Returns
	 *  codes(array)     : (channels,N) counts of close by spots
	 */

	N = coords.length / 4; // number of points
	nc = channels.length;
	codes = newArray(nc * N);

	for (i = 0; i < N; i++) {
	 	ci = coords[4*i];
		xi = coords[4*i+1];
	 	yi = coords[4*i+2];
	 	zi = coords[4*i+3];
	 	for (j = 0; j < N; j++) {
	 		cj =       coords[4*j];
	 		dx = (xi - coords[4*j+1]);
	 		dy = (yi - coords[4*j+2]);
	 		dz = (zi - coords[4*j+3]);
	 		d = dx * dx + dy * dy + dz * dz;
	 		if (d < dmax * dmax) {
	 			codes[cj+i*nc] = codes[cj+i*nc] + 1;
	 		}
	 	}
	 }
	 return codes;
}

function powerSet(n) {
	/* Powerset of n sets as a array (n,2^n)
	 *
	 * List all combination of the first n integers
	 * obtained by using their binary representation.
	 * https://en.wikipedia.org/wiki/Power_set
	 *
	 * parameters
	 *   n (integer): number of items
	 * returns
	 *   set (n,2^n) array with all combinations
	 */
	m = pow(2, n);
	set = newArray(n * m);
	for (i = 0; i < m; i++) {
		bin = toString(IJ.pad(toBinary(i),n));
		for (j = 0; j < n; j++) {
			set[j+n*i] = parseInt(substring(bin, j,j+1));
		}
	}
	return set;
}

function countsCodeBySet(codes, sets, channels) {
	/* Counts the number of events in each combination of channels
	 *
	 * parameters:
	 *  codes (array): (channels, N) counts of close by spots
	 *  sets (array) : (channels, 2^channels) sets
	 *  channels (array): list of channels
	 * returns:
	 *  counts (array): (N,2^channels) array of counts for each sets
	 */
	n = channels.length;
	m = pow(2,n);
	N = codes.length / n;
	counts = newArray(N * m);
	N = codes.length / n;
	for (i = 0; i < N; i++) { // each spots
		for (k = 0; k < m; k++) { // for each set
			acc1 = 0; // sum of neighboors in this combination of channels
			acc2 = 0; // sum of active channels
			test = true; // channel matches the set (channel is >0 when set is 1)
			for (j = 0; j < n; j++) { // for each channel
				acc1 += codes[j + n * i] * sets[j + n * k];
				acc2 += sets[j + n * k];
				if (sets[j + n * k] == 1) {
					test = test && ( codes[j + n * i] > 0);
				}
			}
			if (test) {
				counts[i + N * k] = acc1 / acc2;
			}
		}
	}
	return counts;
}

function aggregateCountsPerSet(counts, set, channels) {
	/* Count the number of non-zeros sets
	 *
	 * parameters:
	 *  counts (array) : (N,2^channels)  array of counts for each sets and spots
	 *  sets (array) : (channels, 2^channels)
	 *  channels (array) : (channels,)
	 * Result:
	 *  agg (array): 2^channels
	 */
	n = channels.length;
	m = pow(2,n);
	N = counts.length / m;
	agg = newArray(m);
	for (k = 0; k < m; k++) { // for each set
		for (i = 0; i < N; i++) {
			if (counts[i + N * k] > 0) {
				agg[k] = agg[k] + 1;
			}
		}
	}
	return agg;
}

function getSetName(set, channels){
	/*get the name of the set of channels as ch1ch2 */
	n = channels.length;
	str = "";
	for (j = 0; j < n; j++) {
		if (set[j] > 0) {
			str += "Ch" + channels[j];
		}
	}
	return str;
}

function rankSetByCardinality(sets, channels) {
	/* rank the sets by the number of active channels */
	n = channels.length;
	m = pow(2,n);
	card = newArray(m);
	for (k = 0; k < m; k++) {
		c0 = 0;
		for (j = 0; j < n; j++) {
			card[k] = card[k] + n * sets[j + n * k];
			if (sets[j + n * k] > 0) {
				c0 = j;
			}
		}
		card[k] = card[k] + c0;
	}
	rank = Array.rankPositions(card);
	return rank;
}

function showCodes(tbl, codes, channels) {
	/*
	 * codes (array): (channels, N) counts of close by spots
	 */
	Table.create(tbl);
	n = channels.length;
	N = codes.length / n;
	for (i = 0; i < N; i++) {
		for (j = 0; j < n; j++) {
			Table.set(channels[j], i, codes[j + n * i]);
		}
	}
	Table.update;
}

function showCounts(tbl, counts, sets, channels) {
	/* Show counts for each set
	 *
	 * tbl: string
	 * 	table name
	 * counts (array)
	 *    (N,2^channels)  array of counts for each sets and spots
	 * sets (array):
	 *   (channels, 2^channels)
	 * channels (array)
	 */
	Table.create(tbl);
	n = channels.length;
	m = pow(2,n);
	N = counts.length / m;
	rank = rankSetByCardinality(sets, channels);
	for (k = 1; k < m; k++) { // for each set
		ko = rank[k];
		set = Array.slice(sets, n * ko, n * (ko + 1));
		col = getSetName(set, channels);
		for (i = 0; i < N; i++) {
			Table.set(col, i, counts[i + N * ko]);
		}
	}
	Table.update;
}

function appendRecord(agg, sets, channels, feature, specificity, spot_size) {
	/* append a record in the table with results */
	tbl = "MultiChannelsColocalisation.csv";
	if (!isOpen(tbl)) {Table.create(tbl);}
	selectWindow(tbl);
	row = Table.size;
	name = File.getNameWithoutExtension(path);
	n = channels.length;
	m = pow(2,n);
	Table.set("Name", row, path);
	rank = rankSetByCardinality(sets, channels);

	total = 0;
	for (k = 1; k < m; k++) {
		ko = rank[k];
		set = Array.slice(sets, n * ko, n * (ko + 1));
		col = getSetName(set, channels);
		Table.set(col, row, round(agg[ko]));
		total = total + agg[ko];
	}
	Table.set("Number", row, round(total));
	Table.set("Feature", row, feature);
	for (i = 0; i < channels.length; i++) {
		Table.set("size Ch" + channels[i], row, spot_size[i]);
		Table.set("specificity Ch" + channels[i], row, specificity[i]);
	}
	Table.update;
}

function setColorFromLut(channel) {
	Stack.setChannel(channel);
	getLut(reds, greens, blues);
	setColor(0.2*255+0.8*reds[255],0.2*255+0.8*greens[255],0.2*255+0.8*blues[255]);
}

function addOverlay(mode, coords, channels, spot_size) {
	if (mode!=3 && !closeonexit) {
		Stack.getDimensions(dummy, dummy, dummy, slices, dummy);
		setBatchMode("hide");
		Overlay.remove();
		getVoxelSize(dx, dy, dz, unit);
		N = coords.length / 4;
		for (i = 0; i < N; i++) {
			c = coords[4*i];
			x = coords[4*i+1] / dx;
			y = coords[4*i+2] / dy;
			z = round(coords[4*i+3] / dz);
			R = 2 * spot_size[c] / dx;
			setColorFromLut(channels[c]);
			if (slices > 1) {
				Stack.setSlice(z);
			}
			Overlay.drawEllipse(x-R, y-R, 2*R+1, 2*R+1);
			if (slices > 1) {
				Overlay.setPosition(channels[c], z, 0);
			} else {
				Overlay.setPosition(channels[c], 1, 0);
			}
			Overlay.add;
		}
		Overlay.show();
		setBatchMode("exit and display");
	}
}

function exportCoordsToCSV(mode, coords, channels, spot_size) {
	if (savecoordinates) {
		opath = File.getDirectory(path) + File.separator + File.getNameWithoutExtension(path) + "-points.csv";
		coords2Table("points.csv", coords, channels, spot_size);
		selectWindow("points.csv");
		Table.save(opath);
	} else {
		if (mode !=3 && !closeonexit) {
			coords2Table("points.csv", coords, channels, spot_size);
		}
	}
}

function zProjectAndShowROIs(mode, coords, channels, spot_size) {
	/* zproject and show roi for a quick evaluation	*/
	if (dozproject && mode!=3 && !closeonexit) {
		run("Select None");
		Stack.getDimensions(_, _, _, slices, _);
		if (slices > 1) {
			run("Z Project...", "projection=[Max Intensity]");
		}
		addOverlay(mode, coords, channels, spot_size);
	}
}

function main() {
	/* Entry point */
	start_time = getTime();
	mode = getMode();
	channels = parseCSVInt(channels_str);
	specificity = parseCSVFloat(specificity_str);
	spot_size = parseCSVFloat(spot_size_str);
	mask_channel = parseMaskStr(mask_str);
	coords = loadData(mode, channels, feature, specificity, spot_size, mask_channel);
	codes = computeCodes(coords, channels, dmax);
	showCodes("codes.csv", codes, channels);
	sets = powerSet(channels.length);
	counts = countsCodeBySet(codes, sets, channels);
	showCounts("counts.csv", counts, sets, channels);
	agg = aggregateCountsPerSet(counts, sets, channels);
	agg1 = correctAggCounts(agg, sets, channels);
	agg2 = divideBySetCardinality(agg1,sets,channels);
	appendRecord(agg2, sets, channels, feature, specificity, spot_size);
	addOverlay(mode, coords, channels, spot_size);
	exportCoordsToCSV(mode, coords, channels, spot_size);
	zProjectAndShowROIs(mode, coords, channels, spot_size);
	end_time = getTime();
	run("Collect Garbage");call("java.lang.System.gc");
	run("Collect Garbage");call("java.lang.System.gc");
	print("Finished in " + (end_time-start_time)/1000 + " seconds.");
	if (closeonexit) {run("Close All");}
}

main();