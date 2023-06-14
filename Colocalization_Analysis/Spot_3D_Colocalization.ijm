// @File(label="Input",description="Use image to run on current image or test to run on a generated test image",value="image") path
// @Integer(label="Channel A", value=1) channela
// @Double(label="Spot size A [um]", value=0.2) spot_sizea
// @Double(label="Specificity A [log10]", value=1, style=slider, min=0, max=6) pfaa
// @Integer(label="Channel B", value=1) channelb
// @Double(label="Spot size B [um]", value=0.2) spot_sizeb
// @Double(label="Specificity B [log10]", value=1, style=slider, min=0, max=6) pfab
// @Double(label="Proximity threshold [um]",value=0.5) dmax

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
tbl = "Spot3DColocalisation.csv";
mode=0;
if (matches(path,".*test")) {
	generateTestImage();
	channela = 1;
	channelb = 2;
	spot_size = 0.2;
	mode=1;
} else if (matches(path,".*image")) {
	print("[Using opened image]");
	mode=2;
} else {
	run("Close All");
	print("Loading image");
	print(path);
	open(path);
	mode=3;
}

// detect spots
A = detect3DSpots(channela, spot_sizea, pfaa);
B = detect3DSpots(channelb, spot_sizeb, pfab);
NA = A.length/3;
NB = B.length/3;

// find pairs of spots
NAB = colocalization(A, B, dmax);
NBA = colocalization(B, A, dmax);

// record results in table
recordInTable(tbl);

// add overlays
if (mode==1 || mode==2) {
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
	Stack.setActiveChannels(active_channels);
	Stack.setDisplayMode("composite");
} 


if (mode==3) {
	close();
}

end_time = getTime();

print("Finished in " + (end_time-start_time)/1000 + " seconds.");

function setColorFromLut(channel) {
	Stack.setChannel(channel);
	getLut(reds, greens, blues);
	setColor(reds[255],greens[255],blues[255]);		
}

function recordInTable(tbl) {
	if (!isOpen(tbl)) {
	Table.create(tbl);	
	}
	
	selectWindow(tbl);
	row = Table.size;
	Table.set("Name", row, getTitle());
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

function coords2Overlay(coords,r,channel) {
	getVoxelSize(dx, dy, dz, unit);	
	R = r / dx;
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

function detect3DSpots(channel, size, pfa) {
	/*
	 * detect spots and return coordinates in physical units as an array
	 */
	setBatchMode(true);
	run("Select None");
	getVoxelSize(dx, dy, dz, unit);	

	sigma1 = size;	
	sigma2 = 3 * size;
	
	id0 = getImageID;
	
	// compute a difference of Gaussian
	run("Duplicate...", "title=id1 duplicate channels="+channel);	
	run("32-bit");
	id1 = getImageID;	
	//run("Square Root","stack");
	run("Gaussian Blur 3D...", "x="+sigma1/dx+" y="+sigma1/dy+" z="+3*sigma1/dz);
	
	run("Duplicate...", "title=id2 duplicate");
	id2 = getImageID;
	run("Gaussian Blur 3D...", "x="+sigma2/dx+" y="+sigma2/dy+" z="+3*sigma2/dz);	
	imageCalculator("Subtract 32-bit stack",id1,id2);
	selectImage(id2);close();
	
	// local maxima
	selectImage(id1);
	run("Duplicate...", "title=id3 duplicate");
	id3 = getImageID;	
	//run("Maximum 3D...", "x=1"+p+" y="+p+" z="+3*p);
	run("Maximum 3D...", "x=2 y=2 z=2");	
	imageCalculator("Subtract 32-bit stack",id3,id1);
	setThreshold(-0.001, 0.0001);
	run("Make Binary", "method=Default background=Dark black");	
	run("Divide...", "value=255 stack");
	
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
	setThreshold(threshold, max);	
	run("Make Binary", "method=Default background=Dark black");	
	run("Minimum 3D...", "x=1 y=1 z=1");	
	run("Maximum 3D...", "x=1 y=1 z=1");		
	run("Divide...", "value=255 stack");	
	
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
	close();
	
	selectImage(id0);
	Stack.setChannel(channel);
	delta = newArray(dx,dy,dz);	
	N = coords.length / 3;	
	for (j = 0; j < N; j++) {				
		p = blockCoordinateGaussianFit(Array.slice(coords,3*j,3*j+3));
		for (k = 0; k < 3; k++) {
			coords[3 * j + k] = p[k] * delta[k];
		}		
	}
	
	setBatchMode(false);
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
	Stack.getDimensions(width, height, channels, slices, frames);
	n = 5;	
	x = p[0];
	y = p[1];
	z = p[2];	
	for (iter = 0; iter < 2; iter++) {
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
		if (random > 0.1) {		
			Stack.setPosition(1, z, 1);
			setPixel(x,y,5000);
			N1++;
			a_on = true;
		}
		if (random > 0.1) {		
			Stack.setPosition(2, z, 1);
			setPixel(x,y,5000);
			N2++;
			b_on = true;
		} 
		if (a_on && b_on) {
			N12++;
		}		
	}
	run("Maximum 3D...", "x=2 y=2 z=2");
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
