// @File(label="Input",description="Use image to run on current image or test to run on a generated test image",value="image") path
// @String(label="Channels", value="1,2,3") channels_str
// @String(label="Spot Size", value="0.2,0.2,0.2") spot_size_str
// @String(label="Specificity[log10]", value="1,1,1") specificity_str
// @Double(label="Proximity threshold [um]",value=0.5) dmax
// @Boolean(label="Subpixel localization", value=true) subpixel
// @Boolean(label="Close on exit", value=false) closeonexit
// @Boolean(label="Save coordinates", value=false) savecoordinates

/*
 * Detect 3D spots and count the number of pairs
 *
 *
 * Input:
 * - Channel: list of channels
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
 * Jerome for Anna2023
 *
 */

start_time = getTime();

mode = 0; // 1:simulation, 2:active window, 3:coordinates, 4:open image
id = 1;

if (matches(path,".*test")) {
	print("Test image");
	generateTestImage();
	channels_str = "1,2,3,4";
	spot_size_str = "0.2,0.2,0.2,0.2";
	specificity_str = "1,1,1,1";
	id = getImageID();
	basename = "test image";
	mode = 1;
} else if (matches(path,".*image")) {
	print("[Using active image "+getTitle()+"]");
	basename = File.getNameWithoutExtension(getTitle);
	id = getImageID();
	mode = 2;
} else if (matches(path,".*csv")) {
	basename = File.getNameWithoutExtension(path);
	print("Loading csv " + path + " table");
	open(path);
	mode = 3;
} else {
	basename = File.getNameWithoutExtension(path);
	print("Opening image");
	print(path);
	run("Bio-Formats Importer", "open=["+path+"] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
	id = getImageID();
	mode = 4;
}

tbl_col = "MultiSpot3DColocalisation.csv";
tbl_pts = basename + "-points.csv";
channels = parseCSVInt(channels_str);
specificity = parseCSVFloat(specificity_str);
spot_size = parseCSVFloat(spot_size_str);

setBatchMode(true);

// Define the name of the table with coordinates
path_pts = File.getDirectory(path) + File.separator + File.getNameWithoutExtension(path) + "-points.csv";

// Load coordinates or detect spots
if (mode == 3) {
	coords = loadCoordsTable(path);
} else {
	print("Detect 3D spots in image " + getTitle);
	coords = newArray(0);
	for (idx = 0; idx < channels.length; idx++) {				
		current_coords = detect3DSpots(channels, spot_size, specificity, idx, subpixel);
		coords2Table(tbl_pts, current_coords, channels, spot_size);
		coords = Array.concat(coords, current_coords);
	}
}

//counts = counts(coords);
//coloc = colocalization(coords, dmax);
//recordInTable(tbl_col, basename, channels, spot_sizes, specificity, dmax, counts, coloc);

/*
 * 
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
	print("Saving file\n" + path_pts);
	selectWindow(tbl_pts);
	Table.save(path_pts);
}

if (closeonexit) {
	selectWindow(tbl_pts);
	run("Close");
	if (!reload) {
		selectImage(id);
		close();
	}
}
*/

end_time = getTime();
run("Collect Garbage");
run("Collect Garbage");
print("Finished in " + (end_time-start_time)/1000 + " seconds.");
//setBatchMode("exit and display");

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


function setColorFromLut(channel) {
	Stack.setChannel(channel);
	getLut(reds, greens, blues);
	setColor(reds[255],greens[255],blues[255]);
}

//function recordInTable(tbl, name, channels, sizes, specificities, dmax, counts, coloc) {
//	
////	if (!isOpen(tbl)) {
////		Table.create(tbl);
////	}
////	
////	selectWindow(tbl);
////	row = Table.size;
////	Table.set("Image name", row, name);
////	Table.set("Proximity threshold [um]", row, dmax);
////	
////	for (i = 0; i < channels.length; i++) {		
////		Table.set("Spot size "+channels[i]+"[um]", row, sizes[i]);
////		Table.set("Specificity "+channels[i]+"[um]", row, specificities[i]);		
////	}
////	
////	for (i = 0; i < nc; i++) {
////		Table.set("N("+channels[i]+")", row, counts[i]);
////	 	if (counts[i] > 0) {
////		 	for (j = 0; j < nc; j++) {
////		 		Table.set("N("+channels[i]+"-"+channels[j]+")", row, coloc[j+i*nc]);		 		
////		 		Table.set("N("+channels[i]+"-"+channels[j]+")/N("+channels[i]+")", row, coloc[j+i*nc]/ counts[i]);		 		
////		 	}
////		}
////	}
//	Table.update();
//}

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
	if (!isOpen(tbl)) {
		Table.create(tbl);
	}
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

function count_channels(coords) {
	// count the number of channels
	N = coords.length / 4;
	nc = 0;
	for (i = 0 ; i < N; i++) {
		if (coords[4*i] > nc) {
			nc = coords[4*i]+1;
		}
	}
	return nc;
}

function counts(coords) {
	// count the number of points in each channel
	N = coords.length / 4; // number of points	
	nc = count_channels(coords);
	counts = newArray(nc); // #Ci
	for (i = 0; i < N; i++) {
		ci = coords[4*i];
		counts[ci] =  counts[ci] + 1;	
	}
	return counts;
}

function create_coloc_codes(coords, dmax) {
	/* colocalization codes */
	N = coords.length / 4; // number of points	 
	nc = count_channels(coords);
	codes = newArray(nc * N);
	for (i = 0; i < N; i++) {
	 	ci = coords[4*i];
		xi = coords[4*i+1];
	 	yi = coords[4*i+2];
	 	zi = coords[4*i+3];	 		 	
	 	counts[ci] = counts[ci] + 1;
	 	for (j = i+1; j < N; j++) if (i!=j) {
	 		cj =       coords[4*j];
	 		dx = (xi - coords[4*j+1]);
	 		dy = (yi - coords[4*j+2]);
	 		dz = (zi - coords[4*j+3]);
	 		d = dx * dx + dy * dy + dz * dz;
	 		if (d < dmax * dmax) {	 			
	 			codes[ci+i*nc] = codes[ci+i*nc] + 1;
	 		}
	 	}	 	
	 }
	 return codes;
}


function codes_counts(codes, channels) {
	nc = channels.length;
	
	
}

function colocalization(coords, dmax) {
	/*
	 * Count number of proximity event
	 *
	 * coords (array)   : coordinate of the set of points (4,n) c,x,y,z
	 * threshold (number) : distance threshold in physical units
	 *
	 * returns return number of close spots for each pairs of channels
	 */
	 
	 N = coords.length / 4; // number of points	 

	 nc = count_channels(coords);
	 	 
	 // allocate arrays
	 coloc = newArray(nc * nc); // #(Ci & Cj)
	 counts = newArray(nc); // #Ci
	 	
	 // count events
	 for (i = 0; i < N; i++) {
	 	ci = coords[4*i];
		xi = coords[4*i+1];
	 	yi = coords[4*i+2];
	 	zi = coords[4*i+3];
	 	dmin = 1e9;
	 	jmin = -1;
	 	counts[ci] = counts[ci] + 1;
	 	for (j = i+1; j < N; j++) if (i!=j) {
	 		cj =       coords[4*j];
	 		dx = (xi - coords[4*j+1]);
	 		dy = (yi - coords[4*j+2]);
	 		dz = (zi - coords[4*j+3]);
	 		d = dx * dx + dy * dy + dz * dz;
	 		if (d < dmax * dmax) {	 			
	 			coloc[cj+ci*nc] = coloc[cj+ci*nc] + 1;
	 		}
	 	}	 	
	 }
	 
	 return coloc;
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
				coords[4*j+3] = slice*dz;
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
	/* generate a test image with 3d spots and c channels */
	run("Close All");
	w = 200;
	h = 200;
	d = 100;
	n = 100;
	c = 4;
	newImage("Test image", "16-bit grayscale-mode", w, h, c, d, 1);
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
