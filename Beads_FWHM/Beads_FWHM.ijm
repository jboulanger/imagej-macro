// @Integer(label="Channel", value=1) channel
// @Double(label="Probablity of false alarm (log10)", value=1, style=slider, min=0, max=5) pfa
// @Double(label="Wavelength [nm]", value=500, style=slider, min=400, max=1000) wavelength_nm
// @Double(label="Numerical Aperture", value=1.25, min=0, max=1.5) numerical_aperture
// @Double(label="Bead size [nm]", value=100, style=slider, min=0.0, max=1000) bead_size_nm
// @String(label="System", value="default") system

run("Close All"); 
if (nImages==0) {
	generateTestImage(wavelength_nm,numerical_aperture,bead_size_nm);
}

// Close the ROI manager if opened
if (isOpen("ROI Manager")) {selectWindow("ROI Manager");run("Close");}

// Compute the expected size of beads to set the spot size for spot detection 
getPixelSize(unit, pw, ph, pd);	
s = 0.51 * wavelength_nm / numerical_aperture;
s = sqrt(s*s + bead_size_nm*bead_size_nm/16);
spot_size = s / 1000 / pw;

// run the spot detection on the specified channel
detect3DSpots(channel, 2*spot_size, pfa);

// add the particles to the ROI manager
run("Analyze Particles...", "add stack");

// close the detection image
close();

// analyse the size of the spots
analyzeSpotsFWHM(system,wavelength_nm,numerical_aperture,bead_size_nm);

function analyzeSpotsFWHM(system,wavelength_nm,numerical_aperture,bead_size_nm) {
	getPixelSize(unit, pw, ph, pd);	
	if (matches(unit,"microns")||matches(unit,"µm")||matches(unit,"micron"))	{
			ps = pw * 1000;			
	} else {
		ps = pw;
	}
	s = 0.51 * wavelength_nm / numerical_aperture;
	s = sqrt(s*s + bead_size_nm*bead_size_nm/16);
	spot_size = s / ps ;
	gamma = 2.4177;// ratio between the FWHM of the Airy function and the std of the Gaussian function		
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
		x0 = (getValue("XM"))/pw-0.5;
		y0 = (getValue("YM"))/ph-0.5;
		Stack.getPosition(channel, z0, frame);
		k = i0 + i;
		Table.set("System", k, system);
		Table.set("Image", k, getTitle());
		Table.set("Channel",k,channel);
		Table.set("X [px]", k, x0);
		Table.set("Y [px]", k, y0);
		Table.set("Z [px]", k, z0);	
		Table.set("Pixel size [nm]", k, ps);

		if (false) {
			profile = circularAverage(x0,y0,10,L);
			radius = getRadiusValues(profile.length);		
			Fit.doFit("Gaussian (no offset)", radius, profile, newArray(1,0,0.21 * wavelength_nm / numerical_aperture * pw));			
			c = Fit.p(2);
		} else {
			vals = circularFit(x0,y0,10,L, 0.21 * wavelength_nm / numerical_aperture * pw);
			Array.getStatistics(vals, min, max, c, stdDev);		
		}		
		c = c * ps;						
		c = sqrt(c * c - bead_size_nm * bead_size_nm / 16);
		Table.set("STD [nm]", k, c);
		Table.set("FWHM [nm]", k, gamma * c);
		Table.set("Abbe [nm]", k, gamma*c/0.514*0.5);		
		Table.set("Sparrow [nm]", k, gamma*c/0.514*0.47);
		Table.set("Rayleigh [nm]", k, gamma*c/0.514*0.61);
		Table.set("Ratio [%]", k, gamma *c / (0.514 * wavelength_nm / numerical_aperture)*100);	
		Table.set("Wavelength [nm]", k, wavelength_nm);
		Table.set("Numerical aperture", k, numerical_aperture);
		Table.set("Bead size [nm] ", k, bead_size_nm);		
	}
}

function circularFit(x0,y0,N,L,s0) {
	getPixelSize(unit, pw, ph, pd);
	makeOval(x0-L+0.5,y0-L+0.5,2*L,2*L);	
	run("Add Selection...");
	vals = newArray(N);
	for (i = 0; i < N; i++) {
		makeLine(x0-L*cos(PI*i/N), y0-L*sin(PI*i/N),x0+L*cos(PI*i/N), y0+L*sin(PI*i/N));		
		run("Add Selection...");
		profile = getProfile();
		radius = getRadiusValues(profile.length);
		Fit.doFit("Gaussian (no offset)", radius, profile, newArray(1,0,s0));		
		vals[i]= Fit.p(2);		
	}
	run("Select None");
	return vals;
}

function getRadiusValues(N) {	
	radius = newArray(N);			
	for (k = 0; k < N; k++) {
		radius[k] = (k-radius.length/2);			
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
	setBatchMode(true);
	run("Select None");
	sigma1 = size;
	sigma2 = 2 * size;
	p = round(2*size+1);
	id0 = getImageID;
	// compute a difference of Gaussian
	run("Duplicate...", "title=id1 duplicate channels="+channel);	
	id1 = getImageID;	
	run("32-bit");
	run("Gaussian Blur 3D...", "x="+sigma1+" y="+sigma1+" z="+sigma1);
	run("Duplicate...", "title=id2 duplicate");
	id2 = getImageID;
	run("Gaussian Blur 3D...", "x="+sigma2+" y="+sigma2+" z="+sigma2);	
	imageCalculator("Subtract 32-bit stack",id1,id2);
	selectImage(id2);close();
	
	// local maxima
	selectImage(id1);
	run("Duplicate...", "title=id3 duplicate");
	id3 = getImageID;		
	run("Maximum 3D...", "x="+p+" y="+p+" z="+p);
	imageCalculator("Subtract 32-bit stack",id3,id1);	
	setThreshold(0,0);
	run("Make Binary", "method=Default background=Dark black");	
	run("Divide...", "value=255 stack");
	
	// threshold
	selectImage(id1);
	if (nSlices>1) {
		Stack.getStatistics(voxelCount, mean, min, max, stdDev);
	} else {
		getStatistics(area, mean, min, max, stdDev, histogram);
	}
	lambda = -2*norminv(pow(10,-pfa));
	//print("pfa="+pfa+" ->" +pow(10,-pfa) + "->"+lambda);
	threshold = mean + lambda * stdDev;
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
	
	// remove edges	
	if (nSlices>1) {
		setSlice(1);
		run("Select All");
		setBackgroundColor(0, 0, 0);
		run("Cut");
		run("Select None");
		setSlice(nSlices);
		run("Select All");
		run("Cut");
		run("Select None");
		}

	// count	
	if (nSlices>1) {
		Stack.getStatistics(voxelCount, mean, min, max, stdDev);
	} else {
		getStatistics(voxelCount, mean, min, max, stdDev, histogram);
	}	
	rename("Spot detection");
	setBatchMode(false);
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
	
	newImage("Stack", "32-bit grayscale-mode", w, h, 1, d, 1);
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
