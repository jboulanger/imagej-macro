#@Integer (label="Channel 1",  value = 2) ch1
#@Integer (label="Channel 2",  value = 3) ch2
#@Float (label="Smoothing in XY [px]", value = 10) sxy
#@Float (label="Smoothing in Z [px]",  value = 10) sz

/*
 * Compute a local correlation coefficient between 2 channels
 * 
 * The Pearson Correlation Coefficient is measured in Gaussian windows
 * defined by the two sigma parameters.
 * 
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
localPearsonCorrelationCoefficient(ch1, ch2, sxy, sz);
setBatchMode(false);

function localPearsonCorrelationCoefficient(ch1, ch2, sxy, sz) {
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
	run("Fire");

	// Clean up
	selectImage(id1); close();
	selectImage(id2); close();
	selectImage(id3); close();
	selectImage(id4); close();
	
	return id5;			
}


function create_test_image() {
	newImage("Test Image", "8-bit composite-mode", 400, 200, 2, 1, 1);
	Stack.setChannel(1);
	run("Macro...", "  code=v=128*(1+sin(2*3.14159*x/32)) slice");
	Stack.setChannel(2);
	run("Macro...", "  code=v=128*(1+sin(3.14159*(2*x/32+x/w))) slice");
}