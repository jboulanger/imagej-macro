/*
 * Digital Phase contrast
 * 
 * Invert the Transport of Intensity Equation to reconstruct the phase from
 * a bright field image stack.
 * 
 * Input: bright field image stack, or stack + time hyperstack
 * Output: 2D or 2D+time brightfield stack
 * 
 * Jerome Boulanger for Anna Edmundson
 */

setBatchMode("hide");
name = File.getNameWithoutExtension(getTitle());
imageDifference();	
solvePoissonFFT();
postprocess();
rename(name+"-DPC");
run("Collect Garbage");
setBatchMode("exit and display");

function postprocess() {
	run("Subtract Background...", "rolling=200 stack");			
	if (nSlices > 1) {		
		Stack.getStatistics(voxelCount, mean, min, max, stdDev);
	} else {
		getStatistics(area, mean, min, max, std, histogram);
	}	
	setMinAndMax(min, max);
	run("16-bit");
}

function solvePoissonFFT() {
	// Solve Laplacian x = Image	
	id = getImageID();	
	run("Duplicate...", "title=flt "); flt = getImageID();
	run("Macro...", "code=[v=-1/(4.001-2*cos(2*(x-w/2+1)*3.14159/w)-2*cos(2*(y-h/2+1)*3.14159/h))]");	
	selectImage(id);
	run("Custom Filter...", "filter=flt process");
	rename("Reconstructed"); rec= getImageID();
	selectImage(flt); close();
	selectImage(rec);		
}

function imageDifference() {
	// take the difference of the two most extreme images	
	Stack.getDimensions(width, height, channels, slices, frames);
	id = getImageID();
	if (frames==1) {
		run("Duplicate...", "duplicate range=1-1");
	} else {
		run("Duplicate...", "title=A duplicate slices=1-1");
	}
	run("32-bit");	
	a = getImageID();
	selectImage(id);
	if (frames==1) {
		run("Duplicate...", "duplicate range="+slices+"-"+slices+"");
	} else {
		run("Duplicate...", "title=Difference duplicate slices="+slices);
	}
	run("32-bit");	
	b = getImageID();		
	imageCalculator("Subtract 32-bit stack", b, a);
	selectImage(a);
	close();
	selectImage(b);	
}

function testSolvePoisson() {
	run("Close All");
	run("Boats");
	run("32-bit");
	run("Smooth");
	run("Convolve...", "text1=[0 1 0\n1 -4 1\n0 1 0\n] normalize");
	solvePoissonFFT();
}
