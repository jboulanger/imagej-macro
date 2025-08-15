#@String(label="Method", choices={"RL","GM"}) method
#@Float (label="Background", value=0) background
#@Float (label="Lateral [um]", value=0.2) lateral
#@Float (label="Axial [um]", value=0.6) axial
#@Float (label="Zoom", value=1) zoom
#@Integer (label="Iterations", value=3) iterations

/*
 * Deblur with RL or GM method
 * 
 * Jerome Boulanger
 */

function forward(sx,sy,sz,zoom) {
	run("Duplicate...", "title=duplicate");			
	run("Gaussian Blur 3D...", "x=" + zoom*sx + " y=" + zoom*sy + " z=" + zoom*sz);
	run("Size...", "width="+width+" height="+height+" depth="+depth+" constrain average interpolation=Bilinear");
	run("Add...", "value=" + background + " stack");
	return getImageID();
}

function backward() {
	run("Add...", "value=" + background + " stack");
}

function deblur(method, background, lateral, axial, zoom, iterations) {
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
	run("Gaussian Blur 3D...", "x=" + sx + " y=" + sy + " z=" + sz);	
	run("Size...", "width="+zoom*width+" height="+zoom*height+" depth="+zoom*depth+" constrain average interpolation=Bilinear");
	run("Min...", "value=0.5");
	u = getImageID();
	for (iter = 0; iter < iterations; iter++) {	
		// Hu
		selectImage(u);
		run("Duplicate...", "duplicate");			
		run("Gaussian Blur 3D...", "x=" + zoom*sx + " y=" + zoom*sy + " z=" + zoom*sz);
		run("Size...", "width="+width+" height="+height+" depth="+depth+" constrain average interpolation=Bilinear");
		run("Add...", "value=" + background + " stack");
		run("Min...", "value=0.5");
		hu = getImageID();
		imageCalculator("Divide create 32-bit stack", f, hu);
		if (method == "RL") {
			run("Size...", "width="+zoom*width+" height="+zoom*height+" depth="+zoom*depth+" constrain average interpolation=Bilinear");
			run("Gaussian Blur 3D...", "x=" + zoom*sx + " y=" + zoom*sy + " z=" + zoom*sz);
		} else {
			run("Size...", "width="+zoom*width+" height="+zoom*height+" depth="+zoom*depth+" constrain average interpolation=Bilinear");
		}
		run("Median 3D...", "x=1 y=1 z=1");
		ratio = getImageID();
		imageCalculator("Multiply 32-bit stack", u, ratio);
		selectImage(hu); close();
		selectImage(ratio); close();
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
	deblur(method, background, lateral, axial, zoom, iterations);
	setBatchMode("exit and display");
}

run("Close All");
main();