/*
 * 3D Debluring with Gold-Meinel algorithm
 * 
 * Jerome Boulanger 2015
 */

macro "Gold_Deblur3D" {
	if (nImages() == 0) {
		run("Confocal Series (2.2MB)");	
	}
	
	Dialog.create("Deblur");
	getPixelSize(unit, pw, ph, pd);	
	if (matches(unit,"microns")||matches(unit,"µm")||matches(unit,"micron"))	{
		ps = pw * 1000;
		zspacing = pd * 1000;
	} else {
		ps = pw;
		zspacing = pd;
	}	
	Dialog.addNumber("Lateral blur (nm):", 100);
	Dialog.addNumber("Axial blur (nm):", 3000);
	Dialog.addNumber("Zoom :", 1);
	Dialog.addNumber("Number of iterations:", 3);
	Dialog.addNumber("Background:", 100);
	Dialog.show();
	lateral = Dialog.getNumber();
	axial = Dialog.getNumber();
	zoom = Dialog.getNumber();
	number_of_iterations = Dialog.getNumber();
	background = Dialog.getNumber();
	gold_deblur(lateral / ps, lateral / ps, axial / pd, zoom, number_of_iterations, background);
}


function gold_deblur(sx, sy, sz, zoom, iter, background) {		
	/* 
	 * 3D debluring using a Gold Meinel algorithm. 
	 * 
	 */
	setBatchMode(true);
	title = getTitle(); // Get windows title
	Stack.getDimensions(width, height, channels, depth, frames);
	id0 = getImageID();
	// Duplicate the image to create a 32bit version of it	
	run("Duplicate...", "duplicate");
	run("32-bit");
	// and substract its minimum - offset to avoid dividing by 0
	//getRawStatistics(nPixels, mean, min);
	run("Subtract...", "value=" + background + " stack");	
	data = getImageID();	
	// Duplicate the data image to create the initialization
	run("Duplicate...", "duplicate");	
	run("Gaussian Blur 3D...", "x=" + sx + " y=" + sy + " z=" + sz);	
	run("Size...", "width="+zoom*width+" height="+zoom*height+" depth="+zoom*depth+" constrain average interpolation=Bilinear");	
	img1 = getImageID(); // ID of the initial estimate
	for (i = 0; i < iter; i++) {		
		// Compute the blurred version of the estimate
		run("Duplicate...", "title=Blurred duplicate");			
		run("Gaussian Blur 3D...", "x=" + zoom*sx + " y=" + zoom*sy + " z=" + zoom*sz);
		run("Size...", "width="+width+" height="+height+" depth="+depth+" constrain average interpolation=Bilinear");
		img2 = getImageID();
		// Divide the input data (id:data) by the blurred estimate (id:img2)
		imageCalculator("Divide create 32-bit stack", data, img2);		
		run("Size...", "width="+zoom*width+" height="+zoom*height+" depth="+zoom*depth+" constrain average interpolation=Bilinear");
		ratio = getImageID();
		// Close the blurred estimate (id:img2)
		selectImage(img2); close();	
		// Multiply the previous estimate (id:img1) by the ration
		imageCalculator("Multiply 32-bit stack", img1, ratio);
		// Close the ratio (id:ratio)
		selectImage(ratio); close();
		// Select the estimate (id:img1)
		selectImage(img1);
		run("Median 3D...", "x=1 y=1 z=1");
	}	
	rename(appendToFilename(title, "_deblurred"));
	for (c = 1; c <= channels; c++) {
		selectImage(id0); Stack.setChannel(c);
		getLut(reds, greens, blues);
		selectImage(img1); Stack.setChannel(c);
		setLut(reds, greens, blues);
		Stack.setSlice(round(zoom*depth/2));		
		run("Enhance Contrast", "saturated=0.05");
	}	
	//Stack.setDisplayMode("composite");
	setBatchMode(false);
	return img1;
}

// append str before the extension of filename
function appendToFilename(filename, str) {
	i = lastIndexOf(filename, '.');
	if (i > 0) {
		name = substring(filename, 0, i);
		ext = substring(filename, i, lengthOf(filename));
		return name + str + ext;
	} else {
		return filename + str
	}
}