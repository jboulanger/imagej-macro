#@String(label="experimental group", description="record an experiement group for book keeping", value="unknown") condition
#@String(label="preprocessing", description="type of object detection", choices={"LoG","Rolling Ball","None"}) preprocess_type
#@Float(label="scale [px]", description="scale of the object to detect", value=2, min=1, max=10, style="slider") scale
#@Float(label="Specificity",description="set the specificity", value=0.25) specificity
#@String(label="area [px]", description="filter object by area, use min-max to set range",value="5-Infinity") areafilter
#@String(label="circularity", description="filter object by circ., use min-max to set range", value="0.5-1") circfilter
#@String(label="overlay colors", description="comma separated list of colors for control overlay", value="blue,red,green,white") colors_str
#@String(label="channel names", description="comma separated list of value of channel name", value="A,B,C") channel_name_str
#@String(label="channels", description="comma separated list of value of channel", value="1,2") channel_str
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
 
print("\\Clear");
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
setBatchMode("hide");
tblname = "coloc.csv";
channels = parseCSVnumbers(channel_str);
channel_names = getChannelNames(channel_name_str);
colors = getColors(colors_str);
img1 = getImageID();
img2 = preprocess(img1, preprocess_type, scale);
parents = getROIList();
print("Number of Parent ROI: " + parents.length);
thresholds = getChannelThresholds(img2, channels, -specificity);
children = getAllchildrenROI(img2, parents, channels, thresholds, areafilter, circfilter, colors);
print("Number of Children ROI: " + children.length);
if (children.length > 1000) {exit;}
children_channel = getROIChannels(children);
measureColocalization(tblname, img1, img2, channels, thresholds, parents, children, children_channel, onoriginal, channel_names);
selectImage(img2); close();
run("Make Composite");
roiManager("show all without labels");
setBatchMode("exit and display");
print("Done");

function getColors(str) {
	colors = split(str,",");
	Stack.getDimensions(width, height, channels, slices, frames);
	while (colors.length < channels) {
		colors = Array.concat(colors, "yellow");
	}
	return colors;
}

function getChannelNames(str) {
	// get the tags for each channel
	a = split(str,",");
	Stack.getDimensions(width, height, channels, slices, frames);
	while (a.length < channels) {
		a = Array.concat(a, newArray("ch"+(1+tags.length)));
	}
	return a;
}

function parseCSVnumbers(str) {
	// return array of numbers from comma separated string
	a = split(str,",");
	b = newArray(a.length);
	for (i=0;i<a.length;i++){b[i] = parseFloat(a[i]);}
	return b;
}


function closeThis(t){
	/* Close the window with title t if opened */	
	if (isOpen(t)) {
		selectWindow(t);
		run("Close");
	}
}

function generate_test_image() {
	w = 1200; // width of the image
	h = 200; // height of the image
	d = 10; // size of the spots
	N = 5; // number of of roi
	n = 5; // number of spots / roi	
	newImage("HyperStack", "32-bit composite-mode", w, h, 2, 1, 1);
	//rho = newArray(0,0.125,0.25,0.5,0.75,1); // amount of colocalization for each ROI
	rho = Array.getSequence(N);
	setColor(255);
	for (roi = 0; roi < rho.length;roi++) {
		rho[roi] = rho[roi] / (rho.length-1);
		makeRectangle(roi*(w/rho.length)+d, d, (w/rho.length)-2*d, h-2*d);
		Roi.setPosition(1,1,1);
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
	run("Stack to Hyperstack...", "order=xyczt(default) channels=2 slices=1 frames=1 display=Color");
	Stack.setChannel(1); resetMinAndMax();
	Stack.setChannel(2); resetMinAndMax();	
}

function preprocess(id0, preprocess_type, scale) {
	selectImage(id0);	
	if (matches(preprocess_type, "LoG")) {		
		run("Select None");
		run("Duplicate...","title=tmp duplicate");
		run("32-bit");		
		run("Square Root", "stack");
		run("Gaussian Blur...", "sigma="+scale+" stack"); 
		run("Convolve...", "text1=[-1 -1 -1\n-1 8 -1\n-1 -1 -1\n] normalize stack");
	} else if(matches(preprocess_type, "Rolling Ball")) {
		run("Select None");
		run("Duplicate...","title=tmp duplicate");
		run("32-bit");		
		run("Gaussian Blur...", "sigma=0.25 scaled stack");
		run("Subtract Background...", "rolling="+5*scale+" stack");
		run("Median...", "radius="+scale+" stack");
	} else {
		run("Select None");
		run("Duplicate...", "duplicate");		
	}
	run("Stack to Hyperstack...", "order=xyczt(default) channels=2 slices=1 frames=1 display=Color");	
	return getImageID();
}

function getROIList() {
	/* Get the list of index of the ROI starting from 0.*/
	roi = newArray(roiManager("count"));
	roi = Array.getSequence(roiManager("count"));	
	return roi;
}

function getChannelThresholds(img, channels, logpfa) {
	/* Return an array of threshold for each channel
	 * img     : id of the window to process
	 * channel : list of channels (array)
	 * logpfa  : log probability of false alarm
	 */
	selectImage(img);		
	run("Select None"); 
	thresholds = newArray(channels.length);
	for (i = 0; i < channels.length; i++) {
		c = channels[i];		
		Stack.setChannel(c);				
		thresholds[c-1] = defineThreshold(logpfa);		
	}
	return thresholds;
}

function getAllchildrenROI(img, parents, channels, thresholds, areafilter, circfilter, colors) {
	/* Create ROI given the list of thresholds for each channel and return the 
	 * list of indices of ROI belonging to all ROI
	 */
	selectImage(img);
	all_children = newArray();	
	for (i = 0; i < channels.length; i++) {		
		children = addChildROI(parents, channels[i], thresholds[i], areafilter, circfilter, colors[channels[i]-1]);
		all_children = Array.concat(all_children, children);
	}
	roiManager("Show All without labels");
	return all_children;
}

function getROIChannels(rois) {
	/* returns the list of channel of the ROIs from their positions
	 * Input
	 * rois : array of ROI indices
	 * Output: array of channel for each roi
	 */
	roi_channels = newArray(rois.length);
	for (i = 0; i < rois.length; i++) {
		roiManager("select", rois[i]);
		//Roi.getPosition(channel, slice, frame);
		Stack.getPosition(channel, slice, frame);		
		roi_channels[i] = channel;
	}	
	return roi_channels;
}

function addChildROI(parents, channel, threshold, area, circ,color) {
	/* Add regions above threshold to the ROI manager and set the channel
	 * Input:
	 * parents : array of parent roi
	 * area : area filter
	 * circ : circularity filter
	 * channel : channel 
	 */
	run("Select None");
	Stack.getDisplayMode(display_mode);
	Stack.setDisplayMode("grayscale"); // to enable threshold on channel	
	// Select all parents to restrict analuse particle to the part of the image in parents ROI	
	if (parents.length > 1) {
		roiManager("select", parents);
		roiManager("combine");
	} else {
		roiManager("select", parents[0]);
	}
	
	Stack.setChannel(channel);
	getStatistics(area, mean, min, max, std, histogram);
	setThreshold(threshold, max);	
	// add new roi to the ROI manager
	n0 = roiManager("count");
	//run("Analyze Particles...", "size="+area+" pixel circularity="+circ+" add");
	run("Analyze Particles...", "add");
	n1 = roiManager("count");	
	roiManager("show all");
	children = newArray(n1 - n0);
	for (i = 0; i < children.length; i++) {
		roiManager("select", n0+i);
		Roi.setPosition(channel, 1, 1);		
		Roi.setStrokeColor(color);
		roiManager("update");
		children[i] = n0 + i;
	}	
	Stack.setDisplayMode(display_mode);
	return children;
}

function hasROIAny(r1, r2) {
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

function hasROIAll(r1, r2) {
	// Test if ROI r2 belong to ROI r2
	ret = true;
	roiManager("select", r2);
	Roi.getCoordinates(x, y);
	roiManager("select", r1);	
	for (i = 0; i < x.length; i++) {
		if (!Roi.contains(x[i], y[i])) {
			ret = false;
			break;
		}
	}
	return ret;
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

function getchildrenUnionInParents(channels, parents, children, child_channels,union_mask) {	
	/* Union of children in given parent and in either channels
	 * Output:
	 * union : the list of indices in ROI manager of selected roi
	 * union_mask : array where the indices of selected children is 1
	 */
	union = newArray(children.length);
	l = 0;
	for (i = 0; i < children.length; i++) {
		// if they belong one of the two channels and roi k
		test = false;
		for (j = 0; j < channels.length; j++) {
			test = test || (child_channels[i] == channels[j]);
		}
		if (test) {
			if (hasROIAny(parents, children[i])) {	
				union[l] = children[i];
				union_mask[i] = true;
				l++;
			}
		} else {
			union_mask[i] = false;
		}
	}
	union = Array.trim(union,l);
	return union;
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

function objectColoc(ch1, ch2, rois, roi_channels, union_mask) {	
	/* Object based colocalization
	 *  
	 */
	if (rois.length != roi_channels.length) {
		print("In objectColoc rois ("+rois.length+") and roi_channels ("+roi_channels.length+") non matching length.")
	}
	if (rois.length != union_mask.length) {
		print("In objectColoc rois ("+rois.length+") and union_mask ("+union_mask.length+") non matching length.")
	}
	print("ch1:" + ch1 + ", ch2:" + ch2);
	Array.print(roi_channels);
	Array.print(union_mask);
	// Compute a subset for rois and roi_channels
	n = 0;
	for (i = 0; i < union_mask.length; i++) if (union_mask[i]) {
		n++;
	}
	subset = newArray(n);
	channels = newArray(n);
	j = 0;
	for (i = 0; i < union_mask.length; i++) if (union_mask[i]) {
		subset[j] = rois[i];
		channels[j] = roi_channels[i];
		j++;
	}
	print("Subset size " + subset.length );
	N1 = 0;
	N2 = 0;
	N12 = 0;	
	for (i = 0; i < subset.length; i++) {		
		ci = channels[i];
		// count the number of roi in each channel
		if (ci == ch1) {
			N1++;
		} else if (ci == ch2) {
			N2++;
		}
		// count the number of overlap
		for (j = i+1; j < subset.length; j++) {
			cj = channels[j];			
			if (ci!=cj && hasROIAny(subset[i], subset[j])) {
				N12++;
			}
		}
	}
	return newArray(N1,N2,N12);
}

function measureColocalization(tblname, img1, img2, channels, thresholds, parents, children, children_channel, on_original, channel_names) {
	/* Measure the colocalization and report by parents
	 * 
	 */
	tblname = "Colocalization.csv";
	if (!isOpen(tblname)) {
		Table.create(tblname);
	}
	
	if (children.length != children_channel.length) {
		print("Children ("+children.length+") and children_channel ("+children_channel.length+") have different size in measureColocalization");
		
	}	
	
	for (parent_id = 0; parent_id < parents.length; parent_id++) {
		
		roiManager("select", parents[parent_id]);
		roi_area = getValue("Area");
		
		for (ch1_id = 0; ch1_id < channels.length; ch1_id++) {
			
			ch1 = channels[ch1_id];
			
			for (ch2_id = ch1_id+1; ch2_id < channels.length; ch2_id++) {								
			
				ch2 = channels[ch2_id];

				union_mask = newArray(children.length);
				union = getchildrenUnionInParents(newArray(ch1,ch2), parents[parent_id], children, children_channel, union_mask);
				
				if (union.length > 0) {
				
					// measure intensity on processed image
					selectImage(img2);
					roiManager("select", union);
					roiManager("combine");
					
					Stack.setChannel(ch1);
					x1 = getIntensities();
					Stack.setChannel(ch2);
					x2 = getIntensities();
					
					// measure intensity on oringal image + min and max to detect saturation
					selectImage(img1);
					roiManager("select", union);
					roiManager("combine");
								
					Stack.setChannel(ch1);
					m1 = getValue("Max");
					y1 = getIntensities();
					Stack.setChannel(ch2);
					m2 = getValue("Max");
					y2 = getIntensities();
																		
					// Compute the colocalization measures
					if (on_original) {
						pcc = pearson(x1, x2);
						rs  = spearman(x1, x2);
						M   = manders(x1, x2,y1, y2, thresholds[ch1_id], thresholds[ch2_id]);
					} else {
						pcc = pearson(y1, y2);
						rs  = spearman(y1, y2);
						M   = manders(x1, x2, x1, x2, thresholds[ch1_id], thresholds[ch2_id]);
					}

					O = objectColoc(ch1, ch2, children, children_channel, union_mask);
								
					// saving results to table
					selectWindow(tblname);
					row = Table.size();
					Table.set("ROI", row, parent_id+1);
					Table.set("Ch1", row, channel_names[ch1-1]);
					Table.set("Ch2", row, channel_names[ch2-1]);
					Table.set("Pearson", row, pcc[0]);
					//Table.set("Spearman", row, rs);
					Table.set("M1", row, M[0]);
					Table.set("M2", row, M[1]);
					Table.set("MOC", row, M[2]);
					Table.set("D", row, (M[5]-M[3]*M[4]/roi_area)/roi_area);

					Table.set("N1", row,  O[0]);
					Table.set("N2", row,  O[1]);
					Table.set("N12", row, O[2]);					
					Table.set("N12/N1", row, O[2] / O[0]);
					Table.set("N12/N2", row, O[2] / O[1]);
									
					Table.set("Area Object 1 [px]", row, M[3]);
					Table.set("Area Object 2 [px]", row, M[4]);
					Table.set("Area Intersection [px]", row, M[5]);
					Table.set("Intersection vs Area Object 1 [px]", row, M[3]/M[5]);
					Table.set("Intersection vs Area Object 2 [px]", row, M[4]/M[5]);
					
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
				}	else {
					print("empty set");
				}	
			}
		}
	}
	Table.update();
}

