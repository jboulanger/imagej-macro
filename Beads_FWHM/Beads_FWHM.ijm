// @Integer(label="Channel", value=1) channel
// @Double(label="Spot size [px]", value=1, style=slider, min=1, max=4.0) spot_size
// @Double(label="Probablity of false alarm (log10)", value=1, style=slider, min=1, max=4.0) pfa
// @Double(label="Wavelength [nm]", value=500, style=slider, min=400, max=1000) wavelength_nm
// @Double(label="Numerical Aperture", value=1.25, min=0, max=1.5) numerical_aperture
// @Double(label="Bead size [nm]", value=100, style=slider, min=0.0, max=1000) bead_size_nm
// @String(label="System", value="default") system

/* 
 *  Beads FWHM
 *  
 *  Compute the full width half max of beads and other resolution criterion based on a Gaussian fit
 *  of the bead.
 *  
 *  Results are saved into a table
 *  
 *  
 *  Jerome Boulanger 2021 for Nick Barry
 */


if (nImages==0) {
	generateTestImage(wavelength_nm,numerical_aperture,bead_size_nm);
}

Overlay.remove();

// Close the ROI manager if opened
if (isOpen("ROI Manager")) {selectWindow("ROI Manager");run("Close");}

setBatchMode(true);

// run the spot detection on the specified channel
Stack.setChannel(channel);
detect3DSpots(channel, spot_size, pfa);

// add the particles to the ROI manager
run("Analyze Particles...", "add stack");

print("Number of detected beads : " + roiManager("count"));

// close the detection image
close();

// analyse the size of the spots
if  (roiManager("count") == 0) {
	print("No detected beads. Aborting.");
} else {
	analyzeSpotsFWHM(system,wavelength_nm,numerical_aperture,bead_size_nm, spot_size);	
}

setBatchMode(false);

function analyzeSpotsFWHM(system,wavelength_nm,numerical_aperture,bead_size_nm,spot_size) {
	getPixelSize(unit, pw, ph, pd);	
	Stack.getDimensions(width, height, channels, slices, frames);
	if (matches(unit,"microns")||matches(unit,"µm")||matches(unit,"micron"))	{
		ps = pw * 1000;
		zspacing = pd * 1000;
	} else {
		ps = pw;
		zspacing = pd;
	}
	gamma = 2.4177;// ratio between the FWHM of the Airy function and the std of the Gaussian function		
	//gamma = 2.355;// ratio between the FWHM of the Gaussian function and its std
	L = 3*spot_size;
	n = roiManager("count");
	i0 = nResults;	
	if (!isOpen("Beads_FWHM.csv")) {
		Table.create("Beads_FWHM.csv");
		i0 = 0;
	} else {
		selectWindow("Beads_FWHM.csv");
		i0 = Table.size;
	}
	for (i = 0; i < n; i++) {
    	roiManager("select", i);
    	List.setMeasurements(); 	
		x0 = (List.getValue("XM"))/pw-0.5;
		y0 = (List.getValue("YM"))/ph-0.5;
		// Roi has wrong channel info because it was create on a duplicated channel
		Stack.getPosition(_channel, z0, frame);
		Stack.setPosition(channel,z0,frame);
		k = i0 + i;
		Table.set("System", k, system);
		Table.set("Image", k, getTitle());
		Table.set("Channel",k,channel);
		Table.set("Wavelength [nm]", k, wavelength_nm);
		Table.set("Numerical aperture", k, numerical_aperture);
		Table.set("Bead size [nm] ", k, bead_size_nm);
		
		Table.set("X [px]", k, x0);
		Table.set("Y [px]", k, y0);
		Table.set("Z [px]", k, z0);	
		Table.set("Pixel size [nm]", k, ps);
		if (false) {
			profile = circularAverage(x0,y0,10,L);
			radius = getRadiusValues(profile.length,ps);	
			Fit.doFit("Gaussian (no offset)", radius, profile, newArray(1,0,0.21 * wavelength_nm / numerical_aperture * pw));			
			c = Fit.p(2);
		} else {
			vals = circularFit(x0,y0,20,L, 0.21 * wavelength_nm / numerical_aperture,ps);
			Array.getStatistics(vals, min, max, c, stdDev);
			Array.getStatistics(vals, cmin, cmax, c, stdDev);
		}
				
		c = sqrt(c * c - bead_size_nm * bead_size_nm / 16);
	
		Table.set("STD XY [nm]", k, c);
		if (slices > 1) {
			cz = axialFit(x0,y0,z0,zspacing);
			cz = sqrt(cz * cz - bead_size_nm * bead_size_nm / 16);
			Table.set("STD Z [nm]", k, cz);
		}
		Table.set("FWHM [nm]", k, gamma * c);
		Table.set("Abbe [nm]", k, gamma*c/0.514*0.5);		
		Table.set("Sparrow [nm]", k, gamma*c/0.514*0.47);
		Table.set("Rayleigh [nm]", k, gamma*c/0.514*0.61);
		Table.set("Ratio [%]", k, gamma *c / (0.514 * wavelength_nm / numerical_aperture)*100);	
	
		Table.set("STD XY min [nm]", k, cmin);
		Table.set("STD XY max [nm]", k, cmax);
	}
}

function axialFit(x0,y0,z0,spacing) {
	/* Fit a Gaussian function axially */	 
	
	Stack.getDimensions(width, height, channels, slices, frames);
	X = newArray(slices);
	Y = newArray(slices);
	for (k = 1; k <= slices; k++) {		
		Stack.setSlice(k);
		makeOval(x0-2, y0-2, 4, 4);
		Y[k-1] = getValue("Max");
		X[k-1] = (k - z0) * spacing;
	}
	Array.getStatistics(Y, ymin, ymax);
	for (k = 0; k < Y.length; k++) {
			Y[k] = (Y[k] - ymin) / (ymax-ymin);
	}
	Fit.doFit("Gaussian", X, Y, newArray(0,1,0,600));
	run("Select None");
	return Fit.p(3);
}

function circularFit(x0,y0,N,L,s0,ps) {
	/* Fit a Gaussian function on profiles extracted over a set of N lines of length L
	 *  Returns a vector of the standard deviation of the profiles
	*/	 
	getPixelSize(unit, pw, ph, pd);
	makeOval(x0-L+0.5,y0-L+0.5,2*L,2*L);	
	run("Add Selection...");
	vals = newArray(N);
	for (i = 0; i < N; i++) {
		makeLine(x0-L*cos(PI*i/N), y0-L*sin(PI*i/N),x0+L*cos(PI*i/N), y0+L*sin(PI*i/N));		
		run("Add Selection...");
		profile = getProfile();
		Array.getStatistics(profile, ymin, ymax);
		for (k = 0; k < profile.length; k++) {
			profile[k] = (profile[k] - ymin) / (ymax-ymin);
		}
		radius = getRadiusValues(profile.length,ps);
		Fit.doFit("Gaussian", radius, profile, newArray(0,1,0, s0));
		vals[i]= Fit.p(3);		
	}
	run("Select None");
	return vals;
}

function getRadiusValues(N,ps) {
	/* Get the N values of radius sampled every ps */	 
	radius = newArray(N);			
	for (k = 0; k < N; k++) {
		radius[k] = (k - N / 2) * ps;			
	}
	return radius;
}

function circularAverage(x0,y0,N,L) {
	getPixelSize(unit, pw, ph, pd);
	makeOval(x0-L+0.5,y0-L+0.5,2*L,2*L);	
	run("Add Selection...");
	for (i = 0; i < N; i++) {
		makeLine(x0-L*cos(PI*i/N), y0-L*sin(PI*i/N),x0+L*cos(PI*i/N), y0+L*sin(PI*i/N));		
		//run("Add Selection...");
		profile = getProfile();
		if (i==0) {
			avg = newArray(profile.length);
		}
		for (k = 0; k < profile.length && k < avg.length; k++) {
			avg[k] = avg[k] + profile[k];
		}
	}
	Array.getStatistics(avg, ymin, ymax);
	for (k = 0; k < avg.length; k++) {
		avg[k] = (avg[k] - ymin) / (ymax-ymin);
	}
	run("Select None");
	Overlay.show();
	return avg;
}

function detect3DSpots(channel, size, pfa) {	
	Stack.getDimensions(width, height, channels, slices, frames);
	run("Select None");
	sigma1 = size;
	sigma2 = 2 * size;
	p = round(2*size+1);
	id0 = getImageID;
	// compute a difference of Gaussian
	run("Duplicate...", "title=id1 duplicate channels="+channel);
	id1 = getImageID;	
	run("32-bit");
	if (slices==1) {
		run("Gaussian Blur...", "sigma="+sigma1);
	} else {
		run("Gaussian Blur 3D...", "x="+sigma1+" y="+sigma1+" z="+sigma1);
	}
	run("Duplicate...", "title=id2 duplicate");
	id2 = getImageID;
	if (slices==1) {
		run("Gaussian Blur...", "sigma="+sigma1);
	} else {
		run("Gaussian Blur 3D...", "x="+sigma2+" y="+sigma2+" z="+sigma2);	
	}
	imageCalculator("Subtract 32-bit stack",id1,id2);
	selectImage(id2); close();
	
	// local maxima
	selectImage(id1);
	run("Duplicate...", "title=id3 duplicate");
	id3 = getImageID;
	if (slices==1) {
		run("Maximum...", "radius="+p);
	} else {
		run("Maximum 3D...", "x="+p+" y="+p+" z="+p);
	}
	imageCalculator("Subtract 32-bit stack",id3,id1);	
	run("Macro...", "code=v=(v==0) stack");	
	
	// threshold
	selectImage(id1);
	if (slices == 1) {
		getStatistics(area, mean, min, max, stdDev, histogram);		
	} else {
		Stack.getStatistics(voxelCount, mean, min, max, stdDev);
	}		
	lambda = -2*norminv(pow(10,-pfa));	
	threshold = mean + lambda * stdDev;
	print("pfa="+pfa+" ->" +pow(10,-pfa) + "->"+lambda+" -> threshold:"+threshold);		
	run("Macro...", "code=v=(v>"+threshold+") stack");	
	
	if (slices==1) {
		run("Minimum...", "radius=1");
		run("Maximum...", "radius=1");				
	} else {
		run("Minimum 3D...", "x=1 y=1 z=1");
		run("Maximum 3D...", "x=1 y=1 z=1");		
	}
		
	// combine local max and threshold
	imageCalculator("Multiply stack",id1,id3);	
	selectImage(id3);close();
			
	// remove edges		
	if (slices > 1) {
		Stack.setSlice(1);
		run("Select All");
		setBackgroundColor(0, 0, 0);
		run("Cut");
		run("Select None");
		Stack.setSlice(slices);
		run("Select All");
		run("Cut");
		run("Select None");
	}
	setThreshold(0.5,1);
		
	rename("Spot detection");		
	return id1;
}

function norminv(x) {
	return -sqrt(  log(1/(x*x)) - 0.9*log(log(1/(x*x))) - log(2*3.14159) );
}

function generateTestImage(wavelength_nm,numerical_aperture,bead_size_nm) {
	w = 200;
	h = 200;
	d = 50;
	n = 5;
	pw = 0.065; //pixel size in micron	
	s = 0.212646 * wavelength_nm / numerical_aperture; // standard deviation in nm matching the airy function
	s = sqrt(s*s + bead_size_nm*bead_size_nm/16);
	spot_size = s / 1000 / pw ;
	
	print("Spot std in pixel "+ spot_size);
	print("Spot std in nm :" + s);	
	
	newImage("Test Image", "32-bit grayscale-mode", w, h, 1, d, 1);
	run("Properties...", "channels=1 slices="+d+" frames=1 pixel_width="+pw+" pixel_height="+pw+" voxel_depth=0.3");
	Stack.setUnits("micron", "micron", "micron", "s", "graylevel");	
	
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			x = round(w*(0.1+0.8*i/(n-1))) ;
			y = round(h*(0.1+0.8*j/(n-1)));
			z = (d+1)/2;
			setSlice(z);
			setPixel(x,y,10000);		
		}
	}

	run("Gaussian Blur 3D...", "x="+spot_size+" y="+spot_size+" z="+spot_size);	
	run("Add Specified Noise...", "stack standard=1");
	resetMinAndMax();
	if (isOpen("Results")) { selectWindow("Results"); run("Close");}

}
