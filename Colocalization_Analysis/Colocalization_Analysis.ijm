#@String(label="preprocessing", description="type of object detection", choices={"LoG","Rolling Ball","None"}) preprocess_type
#@Float(label="scale [px]", description="scale of the object to detect", value=2, min=1, max=10, style="slider") scale
#@Float(label="Specificity",description="set the specificity", value=0.25) specificity
#@String(label="area [px]", description="filter object by area, use min-max to set range",value="5-Infinity") areafilter
#@String(label="circularity", description="filter object by circ., use min-max to set range", value="0.5-1") circfilter
#@String(label="overlay colors", description="comma separated list of colors for control overlay", value="blue,red,green,white") colstr
#@String(label="channel names", description="comma separated list of value of channel name", value="A B C") tagstr
#@String(label="experimental condition", description="record an experiement group for book keeping", value="unknown") condition
#@Boolean(label="Perform colocalization?", description="enable/disable the colocalization", value=true) docoloc
#@Boolean(label="Colocalize on original?", description="use the original image of the preprocessed image", value=true) onoriginal
	
/*
 * Colocalization analysis for segmented structures
 * 
 * - Preprocess all channel with a laplacian of Gaussian filter or a background substraction
 * - Compute a threshold using a probabilty of false alarm based on the background of the image
 * - Add the ROI to the ROI manager using aera and circularity filter, this allows
 *   to discard unwanted objects.
 * - Compute colocalization coefficient within the segmented ROIs for each channel pair
 * 
 * 
 * Usage
 * -----
 * 
 * - Open an image
 * - Add roi in the ROI manager (draw roi and press t)
 * - Run the macro and set the parameters
 * 	o Preprocessing: choose the type of preprocessing (Blobs (lapacian of Gaussian), Regions (median and rolling ball), None)
 *  o Scale : scale in pixel used in the preprocessing step
 *  o False Alarm: control of the false alarm to set the threshold level
 *  o Area[px] : filter ROI by size
 *  o circularity: filter ROI by circularity
 *  o overlay color: list of colors for the overlay of each channel
 *  o channel tag: list of tags for each channel (no space please) 
 *  o condition: add an experimental condition column to the table
 *  o Perform colocalization: allow to run the preprocessing without the colocalization step
 *  o Colocalize on original: colocalization statistics will be all computed on the original image otherwise the preprocessed image 
 *    is used allowing to remove varying background for example. 
 * - Check the results for each cell and inspect the overlays
 * 
 * Output
 * ------
 * 
 * Table with rows :
 * - name: image name
 * - ROI : index of the ROI 
 * - ch1 : first channel 
 * - ch2 : second channel
 * - pearson: Pearson's correlation coefficient (PCC) sum ((R-mR) (G-mG)) / sum((R-mR)^2) * sum((G-mG)^2)
 * - spearman: Spearman rank correlation coefficient
 * - M1 : 1st Mander's coefficient : sum of ch1 intensity in overlap / sum of ch1 intensity in roi 1
 * - M2 : 2nd Mander's coefficient : sum of ch2 intensity in overlap / sum of ch2 intensity in roi 2
 * - MOC : Mander's overlap coefficient sum (ch1 ch2) / sum(ch1^2) * sum(ch2^2)
 * - Mean 1: mean of the 1st channel
 * - Stddev 1 : standard deviation of 1st channel
 * - Max 1 : maximum in original image channel 1
 * - Mean 2: mean of the 2st channel
 * - Stddev 2 : standard deviation of 2st channel
 * - Max 2 : maximum in original image channel 1
 * - N1 : number of pixel above threshold in channel 1
 * - N2 : number of pixel above threshold in channel 2
 * - N12 : number of pixel above threshold in both channel 1 and 2
 * - D = (N12 - N1 N2/n)/n with n the number of pixels in the ROI
 * 
 * Reference
 * ---------
 * 
 * https://jcs.biologists.org/content/joces/131/3/jcs211847.full.pdf
 * https://onlinelibrary.wiley.com/doi/pdf/10.1111/j.1365-2818.1993.tb03313.x
 * 
 * Jerome Boulanger for Yohei 2021
 */

logpfa = -specificity;

// Initialization

if (nImages==0) {
	closeThis("ROI Manager");
	generate_test_image();
}

if (roiManager("count")==0) {
	waitForUser("ROI Manager is empty. Please add at least one ROI");
}

if (roiManager("count")==0) {
	exit("ROI manager was empty.");
}

print("Processing image "+getTitle()+" with scale:"+ scale + " logpfa:" + logpfa + " filters size:"+areafilter + ", circ:" + circfilter);
tstart = getTime();
setBatchMode("hide");


// get image dimensions for loop over channels and pixel sizes for roi position
Stack.getDimensions(width, height, channels, slices, frames);
getPixelSize(unit, pixelWidth, pixelHeight);

// get the colors as an array of string
colors = split(colstr," ");
while (colors.length < channels) {
	colors = Array.concat(colors, "yellow");
}

// get the tags for each channel
tags = split(tagstr," ");
while (tags.length < channels) {
	tags = Array.concat(tags, newArray("ch"+(1+tags.length)));
}

// retreive the id of the image and name
id0 = getImageID(); name = getTitle();

// clear overlays
Overlay.remove();

// add initial ROI as white
roi0 = newArray(roiManager("count"));
roi0 = Array.getSequence(roiManager("count"));
colorROI(roi0,"white",true);

// Processing

// preprocess the image
if (matches(preprocess_type, "LoG")) {
	id  = preprocess1(scale);
} else if(matches(preprocess_type, "Rolling Ball")) {
	id  = preprocess2(scale);
} else {
	run("Select None");
	run("Duplicate...", "duplicate");
	id = getImageID();
}

if (onoriginal) {
	idcoloc = id0;
} else {
	idcoloc = id;
}

// retreive ROI for each channels (allows to filter out by shape)
roi_channel = newArray(); // store the channel of the ROI
all_childs = newArray();
thresholds = newArray(channels);
for (c = 1; c <= channels; c++) {
	print(" - processing channel " + c + "/" + channels);
	selectImage(id);
	run("Select None");
	run("Duplicate...", "duplicate channels="+c);
	run("Select None");
	thresholds[c-1] = defineThreshold(logpfa);
	getStatistics(area, mean, min, max, std, histogram);
	setThreshold(thresholds[c-1], max);
	// add all ROI within the ROI list
	childs = addChildROI(roi0, areafilter, circfilter);
	// Close the temporary image
	close();
	// add the overlay
	selectImage(id0); Stack.setChannel(c);
	colorROI(childs, colors[c-1], false);
	setColor(colors[c-1]);
	setFont("Arial", 20, "bold");
	Overlay.drawString(tags[c-1], 5, c*20);
	// store the ROI for this channel
	roi_col = newArray(childs.length);
	for (i = 0; i < roi_col.length; i++) { roi_col[i] = c;}
	roi_channel = Array.concat(roi_channel, roi_col);
	all_childs = Array.concat(all_childs, childs);
}

if (docoloc) {
	print(" - Colocolization");
	// Compute the colocalization coefficients
	tblname = "Colocalization.csv";
	if (!isOpen(tblname)) {
		Table.create(tblname);
	}
	
	for (k = 0; k < roi0.length; k++) {
		showProgress(k, roi0.length);
		roiManager("select",roi0[k]);
		Roi.getBounds(rx, ry, rwidth, rheight);
		run("Select None");
		roi_area = rwidth*rheight;
		for (c1 = 1; c1 <= channels; c1++) {
			for (c2 = c1+1; c2 <= channels; c2++) {
				// compute the union of roi in roi0[k] and channel 1 and 2
				union = newArray(all_childs.length);
				l = 0;
				for (i = 0; i < all_childs.length; i++) {
					// if they belong one of the two channels and roi k
					if (roi_channel[i] == c1 || roi_channel[i] == c2) {
						if (hasROI(roi0[k], all_childs[i])){
							union[l] = all_childs[i];
							l++;
						}
					}
				}
				union = Array.trim(union,l);
				// Measure the instensity in both channels
				selectImage(id);
				roiManager("select", union);
				roiManager("combine");
				Stack.setChannel(c1);
				x1 = getIntensities();
				Stack.setChannel(c2);
				x2 = getIntensities();
	
				selectImage(id0);
				run("Restore Selection");
				Stack.setChannel(c1);
				m1 = getValue("Max");
				y1 = getIntensities();
				Stack.setChannel(c2);
				m2 = getValue("Max");
				y2 = getIntensities();
			
				// Compute the colocalization measures
				if (onoriginal) {
					pcc = pearson(x1, x2);
					rs  = spearman(x1, x2);
					M   = manders(x1,x2,y1,y2,thresholds[c1-1],thresholds[c2-1]);
				} else {
					pcc = pearson(y1, y2);
					rs  = spearman(y1, y2);
					M   = manders(x1,x2,x1,x2,thresholds[c1-1],thresholds[c2-1]);
				}
							
				// saving results
				selectWindow(tblname);
				row = Table.size();
				Table.set("Condition", row, condition);
				Table.set("Image", row, name);
				Table.set("ROI", row, k+1);
				Table.set("Ch1", row, tags[c1-1]);
				Table.set("Ch2", row, tags[c2-1]);
				Table.set("Pearson", row, pcc[0]);
				Table.set("Spearman", row, rs);
				Table.set("M1", row, M[0]);
				Table.set("M2", row, M[1]);
				Table.set("MOC", row, M[2]);
				Table.set("Area Object 1 [px]", row, M[3]);
				Table.set("Area Object 2 [px]", row, M[4]);
				Table.set("Area Intersection [px]", row, M[5]);
				Table.set("Mean Ch1 in Object 1", row, M[6]);
				Table.set("Mean Ch2 in Object 2", row, M[7]);
				Table.set("Mean Ch1 in Intersection", row, M[8]);
				Table.set("Mean Ch2 in Intersection", row, M[9]);
				Table.set("Mean Ch1", row, pcc[1]);
				Table.set("Stddev Ch1 ", row, pcc[2]);
				Table.set("Max Ch1", row, m1);
				Table.set("Mean Ch2", row, pcc[3]);
				Table.set("Stddev Ch2", row, pcc[4]);
				Table.set("Max Ch2", row, m2);
				Table.set("D", row, (M[5]-M[3]*M[4]/roi_area)/roi_area);
				Table.set("preprocessing",row, preprocess_type);
				Table.set("scale",row, scale);
				Table.set("logpfa",row, logpfa);
				Table.set("min area",row, areafilter);
				Table.set("use original",row, onoriginal);
			}
		}
	}
}
// Close the processed image
selectImage(id); close();
clearROI(all_childs);
// exit batch mode
setBatchMode("show");
selectImage(id0);
Stack.setChannel(1);
print("Done");



// Preprocessing functions
function preprocess1(scale) {
	run("Select None");
	run("Duplicate...","title=tmp duplicate");
	id = getImageID();
	run("32-bit");
	run("Square Root", "stack");
	run("Gaussian Blur...", "sigma="+scale+" stack"); 
	run("Convolve...", "text1=[-1 -1 -1\n-1 8 -1\n-1 -1 -1\n] normalize stack");
	return id;
}

function preprocess2(scale) {
	run("Select None");
	run("Duplicate...","title=tmp duplicate");
	run("32-bit");
	id = getImageID();
	run("Gaussian Blur...", "sigma=0.25 scaled stack");
	run("Subtract Background...", "rolling="+5*scale+" stack");
	run("Median...", "radius=2 stack");
	return id;
}

function defineThreshold(logpfa) {
	 // Define a threshold based on a probability of false alarm (<0.5)
	 st = laplaceparam();
	 pfa = pow(10,logpfa);
	 if (pfa > 0.5) {
	 	return st[0] + st[1]*log(2*pfa);
	 } else {
	 	return st[0] - st[1]*log(2*pfa);
	 }
}

function laplaceparam() {
	// compute laplace distribution parameters 
	getHistogram(values, counts, 256);
	n = getWidth() * getHeight();
	c = 0;
	for (i = 0; i < counts.length && c <= n/2; i++) {
		c = c + counts[i];
	}
	m = values[i-1];
	b = 0;
	for (i = 0; i < counts.length; i++){
		b = b + abs(values[i] - m) * counts[i];
	}
	return newArray(m,b/n);
}

function closeThis(t){
	if (isOpen(t)){selectWindow(t);run("Close");}
}

// Colocalization functions

function getIntensities() {
	// Get the intensities inside the current ROI as an array
	Roi.getContainedPoints(xpoints, ypoints);
	values = newArray(xpoints.length);
	for (i = 0; i < xpoints.length; i++) {
		values[i] = getPixel(xpoints[i], ypoints[i]);
	}
	return values;		
}


function pearson(x1, x2) {
	// Peason correlation coefficient
	//
	// return an array  0:pcc,1:m1,2:sigma1,3:mu2,4:sigma2
	Array.getStatistics(x1, min1, max1, mu1, sigma1);
	Array.getStatistics(x2, min1, max1, mu2, sigma2);
	n = x1.length;
	sxy = 0;
	for (i = 0; i < x1.length; i++) {
		sxy += (x1[i]-mu1)*(x2[i]-mu2);
	}
	pcc = sxy / ((n-1) * sigma1 * sigma2);
	return newArray(pcc,mu1,sigma1,mu2,sigma2);
}

function manders(x1,x2,y1,y2,t1,t2) {
	// Manders overlap coeffcicients
	// x1: value of preprocess image in channel 1
	// x2: value of preprocess image in channel 1
	
	// t1 and t2 are threshold values
	// return an array: 0:m1, 1:m2, 2:moc, 3:n1, 4:n2, 5:n12 
	n1 = 0;
	n2 = 0;
	n12 = 0;
	s1 = 0;//sum R
	s2 = 0;//sum G
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

function spearman(x1,x2) {
	// https://en.wikipedia.org/wiki/Spearman%27s_rank_correlation_coefficient
	r1 = Array.rankPositions(x1);
	r2 = Array.rankPositions(x2);
	sd2 = 0;
	n = r1.length;
	for (i = 0; i < n; i++) {
		d = i - r2[r1[i]];
		sd2 = sd2 + d*d;
	}
	rho = 1 - 6 * sd2 / (n*(n*n-1));
	return rho;
}

// ROI related functions

function addChildROI(parents, area, circ) {
	// get all region with in parents to the ROI manager
	n0 = roiManager("count");
	if (parents.length > 1) {
		roiManager("select", parents);
		roiManager("combine");
	} else {
		roiManager("select", parents[0]);
	}
	run("Analyze Particles...", "size="+area+" pixel circularity="+circ+" add");
	n1 = roiManager("count");
	roiManager("show none");
	childs = newArray(n1-n0);
	for (i = 0; i < childs.length; i++) {
		childs[i] = n0 + i;
	}
	return childs;
}

function hasROI(r1, r2) {
	// Test if ROI r2 belong to ROI r2
	ret = false;
	roiManager("select", r2);
	Roi.getCoordinates(x, y);
	roiManager("select", r1);	
	for (i = 0; i < x.length; i++) {
		if (Roi.contains(x[i], y[i])) {
			ret = true;
			break;
		}
	}
	return ret;
}

function colorROI(roi,color,addlabel) {
	// Add the ROI as overlay in color with an optional label 
	Stack.getPosition(channel, slice, frame);
	for (i = 0; i < roi.length; i++) {
		roiManager("select", roi[i]);
		Overlay.addSelection(color);
		Overlay.setPosition(1, slice, frame);
		Overlay.show();
		if (addlabel) {
			setColor("white");
			setFont("Arial", 40);
			Overlay.drawString(""+(i+1), getValue("X")/pixelWidth, getValue("Y")/pixelHeight);
			Overlay.add();
			Overlay.show();
		}
	}
}

function clearROI(roi) {
	// remove ROIs from the roi manager
	tmp  = Array.sort(roi);
	for (i = tmp.length-1; i >= 0; i--) {
		if (tmp[i] < roiManager("count")) {
			roiManager("select", tmp[i]);
			roiManager("delete");
		}
	}
}


function generate_test_image() {
	w = 1200;
	h = 200;
	d = 10;
	n = 30;
	newImage("HyperStack", "32-bit composite-mode", w, h, 2, 1, 1);
	//rho = newArray(0,0.125,0.25,0.5,0.75,1); // amount of colocalization for each ROI
	rho = Array.getSequence(10);
	setColor(255);
	for (roi = 0; roi < rho.length;roi++) {
		rho[roi] = rho[roi] / (rho.length-1);
		makeRectangle(roi*(w/rho.length)+d, d, (w/rho.length)-2*d, h-2*d);
		roiManager("Add");
		//roiManager("select", roi);
		rx = getValue("BX");
		ry = getValue("BY");
		rw = getValue("Width");
		rh = getValue("Height");
		for (i=0;i<n;i++) {
			x = rx + (rw-d) * random;
			y = ry + (rh-d) * random;
			t = random;
			if (t < 0.5+rho[roi]/2) {
				Stack.setChannel(1);
				makeOval(x, y, d, d);
				fill();
			}
			if (t > 0.5-rho[roi]/2) {
				Stack.setChannel(2);
				makeOval(x, y, d, d);
				fill();
			}
		}
	}
	run("Select None");
	run("Gaussian Blur...", "sigma=2");
	run("Add Specified Noise...", "standard=5");
	rename("Test Image");
	roiManager("show all");
	Stack.setChannel(1); resetMinAndMax();
	Stack.setChannel(2); resetMinAndMax();
}
