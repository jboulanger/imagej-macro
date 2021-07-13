#@Integer (label="Channel 1",  value = 2) ch1
#@Integer (label="Channel 2",  value = 3) ch2
#@String (label="Method", choices={"PCC","Overlap"}) method
#@Float (label="Smoothing in XY [px]", value = 10) sxy
#@Float (label="Smoothing in Z [px]",  value = 10) sz
#@String(label="LUT", choices={"Fire","Grays","Ice"}) lut
#@Boolean(label="Add colorbar", value=true) iwantcolorbar
/*
 * Compute a local correlation coefficient between 2 channels
 * 
 * The Pearson Correlation Coefficient is measured in Gaussian windows
 * defined by the two sigma parameters.
 * 
 * Jerome Boulanger 2021 for Aly Maklouf
 */

if (nImages==0) {
	create_test_image();
	ch1 = 1;
	ch2 = 2;
}

Stack.getDimensions(width, height, channels, slices, frames);
setBatchMode(true);
if (matches(method,"PCC")) {
	localPearsonCorrelationCoefficient(ch1, ch2, sxy, sz, lut);
} else {
	localOverlapCoefficient(ch1, ch2, sxy, sz, lut);
}
setBatchMode(false);

function localPearsonCorrelationCoefficient(ch1, ch2, sxy, sz,lut) {
	// Compute a peasronn Correlation Coefficient in a Gaussian window
	// PCC = sum(I1 I2) / sqrt(sum(I1^2)*sum(I2^2))
	id0 = getImageID();
	name = getTitle();
	Stack.getDimensions(width, height, channels, slices, frames);
	run("Duplicate...", "title=id1 duplicate channels="+ch1);
	selectWindow("id1");
	id1 = getImageID();
	run("32-bit");
	m1 = getValue("Mean");
	run("Subtract...", "value="+m1+" stack");

	selectImage(id0);
	run("Duplicate...", "title=id2 duplicate channels="+ch2);
	selectWindow("id2");
	id2 = getImageID();
	run("32-bit");
	m2 = getValue("Mean");
	run("Subtract...", "value="+m2+" stack");
	
	// Multiply the two channel
	imageCalculator("Multiply create 32-bit stack",id1,id2);
	rename("id3");
	id3 = getImageID();

	// Smoothing with a Gaussian kernel
	selectImage(id1); 
	run("Square","stack");
	if (slices > 1) {
		run("Gaussian Blur 3D...", "x="+sxy+" y="+sxy+" z="+sz);	
	} else {
		run("Gaussian Blur...", "sigma="+sxy+" stack");
	}
	selectImage(id2); 
	run("Square","stack");
	if (slices > 1) {
		run("Gaussian Blur 3D...", "x="+sxy+" y="+sxy+" z="+sz);	
	} else {
		run("Gaussian Blur...", "sigma="+sxy+" stack");
	}
	selectImage(id3); 
	if (slices>1) {
		run("Gaussian Blur 3D...", "x="+sxy+" y="+sxy+" z="+sz);	
	} else {
		run("Gaussian Blur...", "sigma="+sxy+" stack");
	}
	
	// Compute denominator (id1 = id1 x id2)
	imageCalculator("Multiply create 32-bit stack",id1,id2);
	rename("id4");
	id4 = getImageID();
	run("Square Root","stack");

	// Compute the ratio sum(id1 id2) / 
	imageCalculator("Divide create 32-bit stack",id3,id4);
	rename("Colocalization of "+name+" ch"+ch1+"-ch"+ch2);
	id5 = getImageID();
	run(lut);
	setMinAndMax(-1,1);
	if (iwantcolorbar) {
		addColorBar(lut,-1,1);
	}
	
	// Clean up
	selectImage(id1); close();
	selectImage(id2); close();
	selectImage(id3); close();
	selectImage(id4); close();

	run("Select None");
	return id5;			
}

function localOverlapCoefficient(ch1, ch2, sxy, sz,lut) {
	// Compute an overlap coefficient in a Gaussian window
	// M12 = sum(I1 x I2) / sum(I1)
	id0 = getImageID();
	name = getTitle();
	Stack.getDimensions(width, height, channels, slices, frames);
	run("Duplicate...", "title=id1 duplicate channels="+ch1);
	selectWindow("id1");
	id1 = getImageID();
	run("32-bit");
	
	selectImage(id0);
	run("Duplicate...", "title=id2 duplicate channels="+ch2);
	selectWindow("id2");
	id2 = getImageID();
	run("32-bit");
	
	// Multiply the two channel
	imageCalculator("Multiply create 32-bit stack",id1,id2);
	rename("id3");
	id3 = getImageID();

	// Smoothing with a Gaussian kernel
	selectImage(id1); 
	if (slices > 1) {
		run("Gaussian Blur 3D...", "x="+sxy+" y="+sxy+" z="+sz);	
	} else {
		run("Gaussian Blur...", "sigma="+sxy+" stack");
	}
	
	selectImage(id3); 
	if (slices>1) {
		run("Gaussian Blur 3D...", "x="+sxy+" y="+sxy+" z="+sz);	
	} else {
		run("Gaussian Blur...", "sigma="+sxy+" stack");
	}
	
	// Compute the ratio sum(id1 id2) /  sum(id1)
	imageCalculator("Divide create 32-bit stack",id3,id1);
	rename("Colocalization of "+name+" ch"+ch1+"-ch"+ch2);
	id4 = getImageID();

	run(lut);
	Stack.getStatistics(voxelCount, mean, min, max, stdDev);
	setMinAndMax(min,max);
	if (iwantcolorbar) {
		addColorBar(lut,min,max);
	}
	
	// Clean up
	selectImage(id1); close();
	selectImage(id2); close();
	selectImage(id3); close();
	run("Select None");
	return id4;			
}



function create_test_image() {
	newImage("Test Image", "8-bit composite-mode", 400, 200, 2, 1, 1);
	Stack.setChannel(1);
	run("Macro...", "  code=v=128*(1+sin(2*3.14159*x/32)) slice");
	Stack.setChannel(2);
	run("Macro...", "  code=v=128*(1+sin(3.14159*(2*x/32+x/w))) slice");
}


function addColorBar(lut,colorbarmin,colorbarmax) {
	// add a colorbar with min/max given by the array values.
	id0=getImageID();
	run("RGB Color");
	W = getWidth;
	H = getHeight;
	w = round(0.025*W);
	h = round(0.7*H);
	fs = 20; // font size
	img = getImageID();	
	newImage("Colorbar", "32-bit white", w, h, 1);
	cb = getImageID();
	run("Macro...", "code=v=255*(h-y)/h");
	run(lut);
	selectImage(img);
	run("Add Image...", "image=Colorbar x="+W-w-3*fs-10+" y="+(H-h)/2+" opacity=90 zero");
	run("Select All");
	nticks=5;
	for (i = 0; i < nticks; i++) {
		setFont("Arial", fs);
		setColor("yellow");
		str = d2s(colorbarmin+i/(nticks-1)*(colorbarmax-colorbarmin), 1);
		Overlay.drawString(str, W-3*fs-5, (H+h)/2-i/(nticks-1)*h);
	}
	Overlay.show();
	selectImage(cb);close();
	selectImage(id0);
}