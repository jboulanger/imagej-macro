#@String(label="Method", choices={"RL","GM"}) method
#@Float (label="Background", value=0) background
#@Float (label="Lateral [um]", value=200) lateral
#@Float (label="Axial [um]", value=600) axial
#@Float (label="Zoom", value=1) zoom
#@Integer (label="Iterations", value=3) iterations
#@Float(label="Smoothness",value=0.01) smoothness

/*
 * Deblur with RL or GM method
 * 
 * Jerome Boulanger
 */

function hessian(s) {
	if (s > 0) {
		id = getImageID();
		run("FeatureJ Hessian", "smallest smoothing=1.0");
		run("Multiply...", "value="+(-s));
		h = getImageID();
		imageCalculator("Add 32-bit stack", id, h);
		selectImage(h);close();
		selectImage(id);
	}
}

function deblur(method, background, lateral, axial, zoom, iterations, smoothness) {
	f = getImageID();
	Stack.getDimensions(width, height, channels, depth, frames);
	getPixelSize(unit, pw, ph, pd);	
	if (matches(unit,"microns")||matches(unit,"µm")||matches(unit,"micron"))	{
		ps = pw * 1000;
		zspacing = pd * 1000;
	} else {
		ps = pw;
		zspacing = pd;
	}
	sx = lateral / ps;
	sy = lateral / ps;
	sz = axial / zspacing;
	print(sx,sy,sz);
	// init u = Htf
	run("Duplicate...", "title="+appendToFilename(getTitle,"-deblured")+" duplicate");
	run("32-bit");
	run("Subtract...", "value=" + background + " stack");
	if (depth > 1) {
		run("Gaussian Blur 3D...", "x=" + sx + " y=" + sy + " z=" + sz);
		if (zoom !=1) {
			run("Size...", "width="+zoom*width+" height="+zoom*height+" depth="+zoom*depth+" constrain average interpolation=Bilinear");
		}
	} else {
		run("Gaussian Blur...", "sigma="+sx);
		if (zoom !=1) {
			run("Size...", "width="+zoom * width+" height="+zoom * height+" constrain average interpolation=Bilinear");
		}
	}
	run("Min...", "value=0.1");
	u = getImageID();
	for (iter = 0; iter < iterations; iter++) {	
		showProgress(iter, iterations);
		// Hu
		selectImage(u);
		run("Duplicate...", "duplicate");
		if (depth > 1) {
			run("Gaussian Blur 3D...", "x=" + zoom*sx + " y=" + zoom*sy + " z=" + zoom*sz);
			if (zoom !=1) {
				run("Size...", "width="+width+" height="+height+" depth="+depth+" constrain average interpolation=Bilinear");
			}
		} else {
			run("Gaussian Blur...", "sigma=" + zoom*sx);
			if (zoom !=1) {
				run("Size...", "width="+width+" height="+height+" constrain average interpolation=Bilinear");
			}
		}
		run("Add...", "value=" + background + " stack");
		run("Min...", "value=0.5");
		hu = getImageID();
		imageCalculator("Divide create 32-bit stack", f, hu);
		if (method == "RL") {
			if (depth > 1) {
				if (zoom !=1) {
					run("Size...", "width="+zoom*width+" height="+zoom*height+" depth="+zoom*depth+" constrain average interpolation=Bilinear");
				}
				run("Gaussian Blur 3D...", "x=" + zoom*sx + " y=" + zoom*sy + " z=" + zoom*sz);
			} else {
				if (zoom !=1) {
					run("Size...", "width="+zoom*width+" height="+zoom*height+" constrain average interpolation=Bilinear");
				}
				run("Gaussian Blur...", "sigma="+zoom*sx);
			}
		} else {
			if (zoom !=1) {
				if (depth > 1) {
					run("Size...", "width="+zoom*width+" height="+zoom*height+" depth="+zoom*depth+" constrain average interpolation=Bilinear");
				} else {
					run("Size...", "width="+zoom*width+" height="+zoom*height+" constrain average interpolation=Bilinear");
				}
			}
		}
		ratio = getImageID();
		imageCalculator("Multiply 32-bit stack", u, ratio);
		selectImage(hu); close();
		selectImage(ratio); close();
		selectImage(u);
		hessian(smoothness);
		//run("Multiply...", "value=0.7");
		run("Min...", "value=0.001");
	}
}

function appendToFilename(filename, str) {
	// append str before the extension of filename
	i = lastIndexOf(filename, '.');
	if (i > 0) {
		name = substring(filename, 0, i);
		ext = substring(filename, i, lengthOf(filename));
		return name + str + ext;
	} else {
		return filename + str
	}
}

function main() {
	// main entry point
	
	// If there is no image create a test image
	if (nImages == 0) {
		run("Confocal Series");
	}
	setBatchMode("hide");
	deblur(method, background, lateral, axial, zoom, iterations, smoothness);
	setBatchMode("exit and display");
}


main();