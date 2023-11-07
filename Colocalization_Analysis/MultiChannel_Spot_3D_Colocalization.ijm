// @File(label="Input",description="Use image to run on current image or test to run on a generated test image",value="image") path
// @String(label="Channels", value="1,2,3",description="comma separated list of channels") channels_str
// @String(label="Spot Size", value="0.2,0.2,0.2",description="comma separated list of size") spot_size_str
// @String(label="Specificity[log10]", value="1,1,1",description="comma separated list of log pfa") specificity_str
// @Double(label="Proximity threshold [um]",value=0.5,description="define the max distance between two nearby spots") dmax
// @Boolean(label="Subpixel localization", value=true,description="enable subpixel localization") subpixel
// @Boolean(label="Close on exit", value=false,description="close after processing") closeonexit
// @Boolean(label="Save coordinates", value=false) savecoordinates

/*
 * Multiple channel spot 3d detection and colocalization
 * 
 * - The input can be an image stack or a csv file with (c,x,y,z) formats.
 *   Type image for processing the current image or test to run a test as input.
 * 
 * - Channels are comma separated channel indices
 * 
 * A proximity event between for channels (ch1ch2ch3) counts as an event in (ch1), (ch2), (ch3), (ch1ch2), (ch2ch3) and (ch1ch3).
 * 
 */

function getMode() {
	if (matches(path,".*test")) {
		path = "Test Image.tif";
		channels_str = "1,2,3,4";
		spot_size_str = "0.2,0.2,0.2,0.2";
		specificity_str = "1,1,1,1";
		return 1;
	} else if (matches(path,".*image")) {
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

function generateTestImage() {
	/* generate a test image with 3d spots and c channels */
	run("Close All");
	w = 200;
	h = 200;
	d = 100;
	n = 100;
	c = 4;
	name = File.getNameWithoutExtension(path);
	newImage(name, "16-bit grayscale-mode", w, h, c, d, 1);
	run("Macro...", "code=100 stack");
	N1  = 0;
	N2  = 0;
	N12 = 0;
	b   = 10;
	counts = newArray(c*c);
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
		for (j1 = 0; j1 < c; j1++) {
			for (j2 = 0; j2 < c; j2++) {
				if (values[j1] == 1 && values[j2] == 1) {
					counts[j1+c*j2] = counts[j1+c*j2] + 1;
				}
			}
		}
	}
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

function normcdf(x) {
	/*
	 * Normal cumulative density function approximation
	 */
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

function detect3DSpots(channels_list, size_list, pfa_list, channel_idx, subpixel) {
	/*
	 * detect spots and return coordinates in physical units as an array
	 *
	 * Parameter
	 *  channels_list (array)  : list of indices of the channel in the image stack
	 *  size_list (array)
	 *  pfa_list (array) : probability of false alarm in -log10
	 *  channel_idx (int) : index of the channel in the array
	 *  subpixel (boolean): use subpixel localization
	 *
	 *  Returns
	 *   an array (4,n) of coordinates (c,x,y,z) in physical units
	 */

	channel = channels_list[channel_idx];
	size = size_list[channel_idx];
	pfa = pfa_list[channel_idx];

	print(" - detect spots in channel " + channel);

	run("Select None");

	getVoxelSize(dx, dy, dz, unit);
	sigma1 = size;
	sigma2 = 3 * size;
	id0 = getImageID;

	// compute a difference of Gaussian
	run("Duplicate...", "title=id1 duplicate channels="+channel);
	run("32-bit");
	id1 = getImageID;
	run("Square Root","stack");
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
	run("Select None");
	if (nSlices>1) {
		Stack.getStatistics(voxelCount, mean, min, max, stdDev);
	} else {
		getStatistics(area, mean, min, max, stdDev, histogram);
	}
	lambda = -2*normppf(pow(10,-pfa));
	threshold = mean + lambda * stdDev;
	//print("mean:"+mean+", std:"+stdDev+", specicity: "+pfa+", lambda: "+lambda+", threshold: "+threshold);
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
	close();

	// subpixel localization and conversion to physical units
	selectImage(id0);
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

function detectSpotsInAllChannels(id, channels, specificity, spot_size) {	
	setBatchMode("hide");
	print("Detect 3D spots in image " + getTitle);
	coords = newArray(0);
	for (idx = 0; idx < channels.length; idx++) {
		current_coords = detect3DSpots(channels, spot_size, specificity, idx, subpixel);
		coords = Array.concat(coords, current_coords);
	}	
	setBatchMode("exit and display");
	return coords;
}


function loadData(mode, channels, specificity, spot_size) {
	if (mode==1) {
		print("Test image");
		generateTestImage();
		id = getImageID();
		basename = "test image";
		coords = detectSpotsInAllChannels(id, channels, specificity, spot_size);
		//if (isOpen(basename+"-points.csv")) {selectWindow(basename+"-points.csv");run("Close");}
		//coords2Table(basename+"-points.csv", coords, channels, spot_size);
	} else if (mode==2) {
		print("[Using active image "+getTitle()+"]");
		basename = File.getNameWithoutExtension(getTitle);
		id = getImageID();
		coords = detectSpotsInAllChannels(id, channels, specificity, spot_size);
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
		coords = detectSpotsInAllChannels(id, channels, specificity, spot_size);
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
	/* powerset of n sets as a array (n,2^n) */
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

function showCodes(codes, channels) {
	/*
	 * codes (array): (channels, N) counts of close by spots
	 */
	tbl = "Codes.csv";
	Table.create(tbl);
	n = channels.length;	
	N = counts.length / n;	
	for (i = 0; i < N; i++) {		
		for (j = 0; j < n; j++) {
			Table.set(col, i, counts[j + n * i]);
		}
	}
	Table.update;
}

function showCounts(counts, sets, channels) {
	/*
	 * counts (array) : (N,2^channels)  array of counts for each sets and spots
	 */
	tbl = "Counts.csv";
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

function appendRecord(agg, sets, channels, specificity, spot_size) {
	/*append a record in the table with results*/
	tbl = "MultiChannelsColocalisation.csv";
	if (!isOpen(tbl)) {Table.create(tbl);}
	selectWindow(tbl);
	row = Table.size;
	name = File.getNameWithoutExtension(path);
	n = channels.length;
	m = pow(2,n);
	Table.set("Name", row, name);	
	rank = rankSetByCardinality(sets, channels);	
	for (k = 1; k < m; k++) {
		ko = rank[k];
		set = Array.slice(sets, n * ko, n * (ko + 1));
		col = getSetName(set, channels);
		Table.set(col, row, agg[ko]);
	}
	for (i = 0; i < channels.length; i++) {
		Table.set("size Ch" + channels[i], row, spot_size[i]);
		Table.set("specificity Ch" + channels[i], row, specificity[i]);
	}
	Table.update;
}

function setColorFromLut(channel) {
	Stack.setChannel(channel);
	getLut(reds, greens, blues);
	setColor(reds[255],greens[255],blues[255]);
}

function addOverlay(mode, coords, channels, spot_size) {
	if (mode!=3 && !closeonexit) {
		setBatchMode("hide");
		getVoxelSize(dx, dy, dz, unit);
		print(dz);
		N = coords.length / 4;
		for (i = 0; i < N; i++) {
			c = coords[4*i];
			x = coords[4*i+1] / dx;
			y = coords[4*i+2] / dy;
			z = round(coords[4*i+3] / dz);			
			R = 2 * spot_size[c] / dx;
			setColorFromLut(channels[c]);
			Stack.setSlice(z);
			Overlay.drawEllipse(x-R, y-R, 2*R+1, 2*R+1);			
			Overlay.setPosition(0, z, 0);			
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

function main() {

	start_time = getTime();

	mode = getMode();
	channels = parseCSVInt(channels_str);
	specificity = parseCSVFloat(specificity_str);
	spot_size = parseCSVFloat(spot_size_str);
	coords = loadData(mode, channels, specificity, spot_size);
	codes = computeCodes(coords, channels, dmax);
	sets = powerSet(channels.length);
	counts = countsCodeBySet(codes, sets, channels);
	//showCounts(counts, sets, channels);
	agg = aggregateCountsPerSet(counts, sets, channels);
	appendRecord(agg, sets, channels, specificity, spot_size);
	addOverlay(mode, coords, channels, spot_size);
	exportCoordsToCSV(mode, coords, channels, spot_size);
	
	end_time = getTime();
	run("Collect Garbage");
	run("Collect Garbage");
	print("Finished in " + (end_time-start_time)/1000 + " seconds.");
	if (closeonexit) {run("Close All");}
}

main();