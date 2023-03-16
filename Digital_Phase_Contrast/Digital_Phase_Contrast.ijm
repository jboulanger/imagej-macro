//@String(label="Difference mode",choices={"Planes","Stack"},style="radioButtonHorizontal") diff_mode
//@Float(label="Regularization",value=1,min=1,max=10,style=slider,description="power of tens") regularization
//@Boolean(label="Post Process", value=true) do_postprocess
/*
 * Digital Phase contrast
 * 
 * Invert the Transport of Intensity Equation to reconstruct the phase from
 * a bright field image stack.
 * 
 * The TIE equation is Delta phi = - (lamda/2pi) I dI/dz
 * 
 * We assume I to be constant
 * 
 * Input: bright field image stack, or stack + time hyperstack
 * Output: 2D or 2D+time brightfield stack
 * 
 * The macro takes the difference of the two most extreme planes of each z-stack
 * and sovle the Poisson equation Delta f = -dI/dz. A post processing applies
 * a rolling ball to remove phase trends in the reconstructed image.
 * 
 * If no image is provided, a test image is generated
 * 
 * Jerome Boulanger for Anna Edmundson
 */
 

if (nImages==0) {	
	testSolvePoisson();
	showMessage("Test 1");
	testImage();
}

setBatchMode("hide");
name = File.getNameWithoutExtension(getTitle());
imageDifference(diff_mode);	
solvePoissonFFT(regularization);
if (do_postprocess) {
	postprocess();
}
rename(name+"-DPC");
run("Collect Garbage");
setBatchMode("exit and display");

function postprocess() {
	run("Subtract Background...", "rolling=200 sliding stack");			
	if (nSlices > 1) {		
		Stack.getStatistics(voxelCount, mean, min, max, stdDev);
	} else {
		getStatistics(area, mean, min, max, std, histogram);
	}	
	setMinAndMax(min, max);
	run("16-bit");
}

function solvePoissonFFT(regularization) {
	// Solve Laplacian x = Image	
	id = getImageID();	
	run("Duplicate...", "title=flt "); flt = getImageID();	
	run("Macro...", "code=[v=-1/(4+"+pow(10,-regularization)+"-2*cos(2*(x-w/2)*3.14159/w)-2*cos(2*(y-h/2)*3.14159/h))]");	
	selectImage(id);
	run("Custom Filter...", "filter=flt process");
	rename("Reconstructed"); rec= getImageID();
	selectImage(flt); close();
	selectImage(rec);		
}

function imageDifference(mode) {
	/* Image difference
	 *  mode : if 1 takes only the two most extreme images, if 0 finite difference in z
	 * Returns either an plane or a stack 
	 */
	Stack.getDimensions(width, height, channels, slices, frames);
	id = getImageID();
	if (matches(mode, "Planes")) {
		if (frames==1) {
			run("Duplicate...", "duplicate range=1-1");
		} else {
			run("Duplicate...", "title=A duplicate slices=1-1");
		}
	} else {
		if (frames==1) {
			run("Duplicate...", "duplicate range=0-" + (slices - 1));
		} else {
			run("Duplicate...", "title=Difference duplicate slices=0" + (slices - 1));
		}
	}
	run("32-bit");	
	a = getImageID();
	selectImage(id);
	if (matches(mode, "Planes")) {
		if (frames==1) {
			run("Duplicate...", "duplicate range="+slices+"-"+slices+"");
		} else {
			run("Duplicate...", "title=Difference duplicate slices="+slices);
		}
	} else {		
		if (frames==1) {
			run("Duplicate...", "duplicate range=2-" + slices);
		} else {
			run("Duplicate...", "title=A duplicate slices=2-" + slices);
		}
	}	
	run("32-bit");	
	b = getImageID();		
	imageCalculator("Subtract 32-bit stack", b, a);	
	selectImage(a);	
	close();
	selectImage(b);	
	return b;
}


function testSolvePoisson() {
	run("Close All");
	run("Boats");
	run("32-bit");
	run("Smooth");
	run("Convolve...", "text1=[0 1 0\n1 -4 1\n0 1 0\n] normalize");
	solvePoissonFFT(regularization);
}


function testImage() {
	// generate a test image
	run("Close All");
	run("Boats");
	run("Grays");
	run("32-bit");
	run("Gaussian Blur...", "sigma=1");	
	id0 = getImageID();
	
	// create a phase image
	makeRectangle(59, 132, 256, 256);
	run("Duplicate...", "title=phase ");
	phi0 = getImageID();
	
	// create an intensity image
	selectImage(id0);
	makeRectangle(355, 128, 256, 256);
	run("Duplicate...", "title=intensity ");
	setColor(1);	fill();
	I0 = getImageID();
	
	// compute the laplacian of the phase
	selectImage(phi0);
	run("Duplicate...", " ");
	run("Convolve...", "text1=[0 1 0\n1 -4 1\n0 1 0\n]");
	deltaphi = getImageID(); rename("delta phase");
	
	//compute an image 
	imageCalculator("Multiply create 32-bit", I0, deltaphi);
	I1 = getImageID();
	run("Multiply...", "value=-1");
	imageCalculator("Add 32-bit", I1, I0);
	rename("intensity plane 0");
	
	imageCalculator("Multiply create 32-bit", I0, deltaphi);
	I2 = getImageID();
	imageCalculator("Add 32-bit", I2, I0);
	rename("intensity plane 1");
	
	run("Images to Stack", "  title=[intensity plane] use");
	stack = getImageID();
	selectImage(deltaphi); close();
	selectImage(id0); close();
	return stack;
}


