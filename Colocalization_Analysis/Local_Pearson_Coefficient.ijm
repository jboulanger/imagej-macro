#@Integer (label="Channel 1",  value = 2) ch1
#@Integer (label="Channel 2",  value = 3) ch2
#@Float (label="Smoothing in XY [px]", value = 10) sxy
#@Float (label="Smoothing in Z [px]",  value = 10) sz

/*
 * Compute a local correlation coefficient between 2 channels
 * 
 * 
 * Jerome Boulanger for Aly Maklouf
 */

Stack.getDimensions(width, height, channels, slices, frames);
setBatchMode(true);
localPearsonCorrelationCoefficient(2, 3, sxy, sz);
setBatchMode(false);

function localPearsonCorrelationCoefficient(ch1, ch2, sxy, sz) {
	// Compute a peasronn Correlation Coefficient in a Gaussian window
	// PCC = sum(I1 I2) / sqrt(sum(I1^2)*sum(I2^2))
	id0 = getImageID();
	ch1 = 2;
	ch2 = 3;
	sxy = 20;
	sz = 2;
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
	run("Gaussian Blur 3D...", "x="+sxy+" y="+sxy+" z="+sz);	
	selectImage(id2); 
	run("Square","stack");
	run("Gaussian Blur 3D...", "x="+sxy+" y="+sxy+" z="+sz);
	selectImage(id3); 
	run("Gaussian Blur 3D...", "x="+sxy+" y="+sxy+" z="+sz);

	// create a mask
	/*
	selectImage(id3);
	run("Duplicate...", "title=mask duplicate");
	selectWindow("mask");
	mask = getImageID();	
	Stack.getStatistics(voxelCount, mean, min, max, stdDev);
	setThreshold(10, max);
	run("Make Binary", "method=Default background=Dark black");
	run("Divide...", "value=255 stack");
	run("32-bit");
	*/
	
	// Compute denominator (id1 = id1 x id2)
	imageCalculator("Multiply create 32-bit stack",id1,id2);
	rename("id4");
	id4 = getImageID();
	run("Square Root","stack");

	// Compute the ratio sum(id1 id2) / 
	imageCalculator("Divide create 32-bit stack",id3,id4);
	rename("Colocalization");
	id5 = getImageID();

	//imageCalculator("Multilpy 32-bit stack",id5,mask);

	// Clean up
	//selectImage(id1); close();
	//selectImage(id2); close();
	//selectImage(id3); close();
	//selectImage(id4); close();
		
}