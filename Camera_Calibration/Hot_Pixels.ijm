/*
 * Correct Hot pixels in a stack
 */
name = appendToFilename(getTitle,"-hpc");
run("Duplicate...", "title=["+name+"] duplicate");	
setBatchMode("hide");
correct_hotpixels();
setBatchMode("show");


function correct_hotpixels() {	
	id = getImageID;	
	run("Select None");
	// create the median image
	selectImage(id);
	run("Duplicate...", "title=med duplicate");	
	med = getImageID();
	//run("Median...", "radius=3 stack");
	run("Median 3D...", "x=3 y=3 z=2");
	// create a mask
	selectImage(id);
	run("Duplicate...", "title=mask duplicate");	
	run("32-bit");
	mask = getImageID();
	run("Convolve...", "text1=[0 -1 0\n-1 4 -1\n0 -1 0]  stack");	
	imageCalculator("Divide stack", mask, med);
	if (nSlices > 1)	{
		Stack.getStatistics(voxelCount, mean, min, max, std);
	} else {
		getStatistics(area, mean, min, max, std, histogram);
	}
	t = 2 * std;	
	run("Macro...", "code=v=abs(v)<"+t+" stack");	
	imageCalculator("Multiply stack", id, mask);
	// invert the mask
	selectImage(mask);
	run("Macro...", "code=v=1-v stack");	
	imageCalculator("Multiply stack", med, mask);
	// add the two images
	imageCalculator("Add stack", id, med);	
	selectImage(mask); close();
	selectImage(med); close();
	selectImage(id);	
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
	
function correct_hotpixels_old(){
	run("32-bit");	
	id = getImageID();
	run("Duplicate...", "duplicate");
	run("32-bit");		
	Stack.getDimensions(width, height, channels, slices, frames);
	run("Stack to Hyperstack...", "order=xyczt(default) channels=1 slices="+channels*slices*frames+" frames=1 display=Color");	
	tmp = getImageID();
	E = getE(tmp);
	V = getV(tmp);	
	selectImage(E);
	offset = 0.5 * getValue("Min");
	run("Subtract...", "value="+offset+" stack");
	imageCalculator("Divide stack", V, E);
	run("Z Project...", "projection=Median");
	G = getImageID();rename("G");
	selectImage(tmp);close();
	selectImage(E);close();
	selectImage(V);close();
	selectImage(id);
	run("Subtract...", "value=100 stack");
	imageCalculator("Divide stack", id, G);
	selectImage(id);
	run("Add...", "value=100 stack");
}

function getE(id) {
	selectImage(id);
	run("Duplicate...", "duplicate");	
	run("Median 3D...", "x=0 y=0 z=5");
	rename("E");
	return getImageID();
}

function getV(id) {
	selectImage(id);
	run("Duplicate...", "title=res duplicate");	
	res = getImageID;
	run("Convolve...", "text1=[0 -0.1826 0\n-0.1826 0.7303 -0.1826\n0 -0.1826 0]  stack");	
	
	run("Duplicate...", "title=med duplicate");
	med = getImageID;	
	run("Median 3D...", "x=0 y=0 z=5");	
	//run("Subtract Background...", "rolling=5");
	
	imageCalculator("Difference stack", res, med);
	selectImage(res);
	run("Abs", "stack");	
	run("Median 3D...", "x=0 y=0 z=5");
	run("Multiply...", "value=0.68");
	//run("Subtract Background...", "rolling=5");
	run("Square","stack");
	rename("V");
	selectImage(med);close();
	selectImage(res);
	return res;
}