//@String(label="channel names") channel_names_str
/*
 * Analysis of colocalization coefficient in a tissue
 * 
 * - Segment nuclei using startdist applied toa DAPI channel
 * - Compute a tesselation of the image using the marker with a watershed from morpholibJ
 * - Convert the labels to ROIs
 * - Compute the colocalization coefficients in the ROIs
 * 
 * Requires CSBDeep, StartDist and MorpholibJ
 * 
 * Jerome Boulanger for Leonor 2023
 */
  
 channel_names = parseCSVString(channel_names_str);
 Overlay.remove;
 run("Select None");
 id = getImageID();
 setBatchMode("hide");
 segment_nuclei(id, 1);
 colocalization(id, channel_names);
 addROIsToOverlay(id);
 setBatchMode("exit and display");
 run("Collect Garbage");	
 run("Select None");
 for (c = 1; c <= nSlices; c++) {
 	Stack.setChannel(c);
 	run("Enhance Contrast", "saturated=0.15");
 }
 Stack.setDisplayMode("composite");
 
 function segment_nuclei(id,channel) { 	
 	print("Segmentation [ mem:" + round(IJ.currentMemory()/1e6)+"MB]");
 	selectImage(id);
 	run("Duplicate...", "title=dapi duplicate channels="+channel);
 	nuc = getImageID();
 	print("Running StarDist [ mem:" + round(IJ.currentMemory()/1e6)+"MB]");
 	run("Square Root");
 	run("Command From Macro", "command=[de.csbdresden.stardist.StarDist2D], args=['input':'dapi','modelChoice':'DSB 2018 (from StarDist 2D paper)', 'normalizeInput':'true', 'percentileBottom':'5.0', 'percentileTop':'95.0', 'probThresh':'0.5', 'nmsThresh':'0.4', 'outputType':'Label Image', 'nTiles':'9', 'excludeBoundary':'0', 'roiPosition':'Automatic', 'verbose':'false', 'showCsbdeepProgress':'false', 'showProbAndDist':'false'], process=[false]"); 	

 	print("Creating  a mask [ mem:" + round(IJ.currentMemory()/1e6)+"MB]");
 	mask = createMask(id);

 	print("Running Watershed [ mem:" + round(IJ.currentMemory()/1e6)+"MB]");
 	selectImage("Label Image"); 
 	run("Marker-controlled Watershed", "input=[Label Image] marker=[Label Image] mask=[Mask] compactness=0 calculate use");
 	run("Collect Garbage");	
 	
 	print("Creating ROIs [ mem:" + round(IJ.currentMemory()/1e6)+"MB]");
 	nlabels = getValue("Max");
 	print("Number of regions " + nlabels);
 	if (isOpen("ROI Manager")) {selectWindow("ROI Manager");run("Close");} 	
 	getPixelSize(unit, pixelWidth, pixelHeight);
	a = getWidth() * getHeight() * pixelWidth * pixelHeight;
 	for (i = 1; i <= nlabels; i++) {
 		setThreshold(i,i); 		
 		run("Create Selection");
 		if (getValue("Area") < a) {
 			roiManager("add");
 		}
 	} 	
 	selectImage(nuc); close();
 	selectImage(mask); close();
 	selectImage("Label Image"); close();
 	selectImage("Label Image-watershed"); close();
 	run("Collect Garbage");	
}

function createMask(id) {
	selectImage(id);
	Stack.setChannel(2);
	Stack.setDisplayMode("color");
	run("Duplicate...", "title=Mask");
	run("Median...", "radius=5");
	percentileThreshold(5);	
	return getImageID();
}

function percentileThreshold(p) {
	getHistogram(values, counts, 1000);
	n = 0;
	getDimensions(width, height, channels, slices, frames);
	for (i = 0; i < counts.length && n < p / 100 * width * height; i++) {
		n += counts[i]; 
	}
	setThreshold(values[i], 65535, "raw");
	run("Convert to Mask");
}

function addROIsToOverlay(id) {
	selectImage(id);
	n = roiManager("count");
	for (i = 0; i < n; i++) {
		roiManager("select", i);
		run("Add Selection...");	
		Overlay.setStrokeColor("white");
	}		
}

function parseCSVString(csv) {		
	str = split(csv,",");
	values = newArray(str.length);
	for (i = 0 ; i < str.length; i++) {
		values[i] = String.trim(str[i]);
	}
	Stack.getDimensions(width, height, channels, slices, frames);
	for (i = str.length; i < channels; i++) {
		values[i] = "ch" + i;
	}
	return values;
}

function colocalization(id, names) {
	print("Colocalization [ mem:" + round(IJ.currentMemory()/1e6)+"MB]");
	selectImage(id);
	Stack.getDimensions(width, height, channels, slices, frames);
	n = roiManager("count");
	for (i = 0; i < n; i++) {
		roiManager("select", i);
		Table.set("ROI",i,i+1);
		Table.set("Area",i,getValue("Area"));
		for (c1 = 1; c1 <= channels; c1++) {
			Stack.setChannel(c1);
			I1 = getIntensities();
			for (c2 = c1+1; c2 <= channels; c2++) {
				Stack.setChannel(c2);
				I2 = getIntensities();
				pcc = pearson(I1, I2);				
				Table.set("PCC ["+names[c1-1]+":"+names[c2-1]+"]", i, pcc[0]);
				
			}
		}
	}	
}

function getIntensities() {
	// Get the intensities inside the current ROI as an array
	Roi.getContainedPoints(xpoints, ypoints);
	values = newArray(xpoints.length);
	for (i = 0; i < xpoints.length; i++) {
		values[i] = getPixel(xpoints[i], ypoints[i]);
	}
	return values;
}

function spearman(x1,x2) {
	// https://en.wikipedia.org/wiki/Spearman%27s_rank_correlation_coefficient
	// The ranking function is not consistent with other statistical software
	r1 = Array.rankPositions(x1);
	r2 = Array.rankPositions(x2);
	pcc = pearson(r1,r2);
	return -pcc[0];
}

function pearson(x1, x2) {
	// Peason correlation coefficient
	// return an array  0:pcc,1:mu1,2:sigma1,3:mu2,4:sigma2
	if (x2.length != x1.length ) {
		print("non matching array length in prearson x1:"+x1.length+", x2:"+x2.length);
	}
	Array.getStatistics(x1, min1, max1, mu1, sigma1);
	Array.getStatistics(x2, min1, max1, mu2, sigma2);
	n = x1.length;
	sxy = 0;
	sxx = 0;
	syy = 0;
	for (i = 0; i < x1.length; i++) {
		sxy += (x1[i]-mu1)*(x2[i]-mu2);
		sxx += (x1[i]-mu1)*(x1[i]-mu1);
		syy += (x2[i]-mu2)*(x2[i]-mu2);
	}
	pcc = sxy / (sqrt(sxx)*sqrt(syy));
	return newArray(pcc,mu1,sigma1,mu2,sigma2);
}

function manders(x1,x2,y1,y2,t1,t2) {
	// Manders overlap coeffcicients
	// x1: value of preprocess image in channel 1
	// x2: value of preprocess image in channel 1
	// t1 and t2 are threshold values
	// return an array: 0:m1, 1:m2, 2:moc, 3:n1, 4:n2, 5:n12

	if (x2.length != x1.length || y1.length != x1.length || y2.length != x1.length) {
		print("non matching array length in manders x1:"+x1.length+", x2:"+x2.length+"y1:"+y1.length+", y2:"+y2.length);
	}

	n1 = 0;
	n2 = 0;
	n12 = 0;
	s1 = 0; // sum R
	s2 = 0; // sum G
	s12 = 0;// sum Rcoloc
	s21 = 0;// sum Gcoloc
	moc12 = 0;
	moc1 = 0;
	moc2 = 0;
	for (i = 0; i < x1.length; i++) {
		if (x1[i] > t1) {n1++; s1 += y1[i]; }
		if (x2[i] > t2) {n2++; s2 += y2[i]; }
		if (x1[i] > t1 && x2[i] > t2) { n12++;  s12 += y1[i]; s21 += y2[i];}
		moc12 += y1[i]*y2[i];
		moc1  += y1[i]*y1[i];
		moc2  += y2[i]*y2[i];
	}
	m1 = s12 / s1;
	m2 = s21 / s2;
	moc = moc12 / sqrt(moc1 * moc2);
	return newArray(m1,m2,moc,n1,n2,n12,s1/n1,s2/n2,s12/n12,s21/n12);
}
