#@Integer(label="Signal channel",value=2) ch_obj
#@Float(label="Threshold",value=128) threshold
/*
 * Quantify distribution of intensity inside a region of interest
 * 
 * Usage:
 * 1. Open an image
 * 2. Create ROI and add them in the ROI manager (press t)
 * 3. Open the macro and press "Run"
 * 4. Set the channels and the threshold
 * 
 * Jerome for Felicity 2021
 */


getVoxelSize(dx, dy, dz, unit);
ustr = " ["+unit+"]";

// Add an overlay to show where the threshold was set
Overlay.remove();
Stack.getDimensions(width, height, channels, slices, frames);
Stack.getStatistics(voxelCount, mean, min, max, stdDev);
for (z = 1; z <= slices; z++) {
	Stack.setPosition(ch_obj, z, 1);	
	setThreshold(threshold, max);	
	run("Create Selection");
	print(getValue("Area"));
	if (getValue("Area")  > 0 && getValue("Area") < width*height*dx*dy/2) {
		run("Add Selection...");
		Overlay.setPosition(ch_obj, z, 1);
	}
}

// 
row = nResults;
n = roiManager("count");
Stack.setChannel(ch_obj);
for (i = 0; i < n; i++) {
	roiManager("select",i);
	xyzv = getROIXYZVData(ch_obj);	
	spread = computeROISpreadXYZ(xyzv);
	eigs = symeigen3x3(spread);
	Array.print(spread);
	entropy = computeROIEntropy(xyzv);	
	volume = computeROITotalVolume(xyzv,threshold);
	setResult("Image", row, getTitle);
	setResult("ROI", row, i+1);
	setResult("Channel", row, ch_obj);
	setResult("Spread in X"+ustr, row, sqrt(spread[0]));
	setResult("Spread in Y"+ustr, row, sqrt(spread[3]));
	setResult("Spread in Z"+ustr, row, sqrt(spread[5]));
	setResult("Spread total"+ustr, row, sqrt(spread[0]+spread[3]+spread[5]));
	setResult("Spread major"+ustr, row, sqrt(eigs[2]));
	setResult("Spread minor"+ustr, row, sqrt(eigs[0]));
	setResult("Entropy", row, entropy);
	setResult("Volume above threshold [voxels]", row, volume);
	setResult("Volume above threshold [" + unit+ "^3]", row, volume*dx*dy*dz);
	row++;
}
Stack.setChannel(ch_obj);

function getROIXYZVData(ch) {
	// get x,y,z,intensity from all pixels in the ROI
	Stack.getDimensions(width, height, channels, slices, frames);	
	getVoxelSize(dx, dy, dz, unit);
	Roi.getContainedPoints(x, y);
	n = 4 * x.length * slices;
	ret = newArray(n);	// x,y,z,intensity
	for (z = 1; z <= slices; z++) {
		Stack.setPosition(ch, z, 1);		
		for (i = 0; i < x.length; i++) {		
			idx = 4 * (i + x.length * (z-1));
			ret[idx]	= x[i]*dx;
			ret[idx+1]	= y[i]*dy;
			ret[idx+2]	= z*dz;
			ret[idx+3]	= getValue(x[i],y[i]);
		}	
	}	
	return ret;
}

function computeROISpreadXYZ(xyzv) {
	/* Compute the spread in x,y and z for the xyzv input	
	 *  
	 */
	swx = 0;
	swxx = 0;
	swy = 0;
	swyy = 0;
	swz = 0;
	swzz = 0;
	sw = 0;	
	swxy = 0;
	swxz = 0;
	swyz = 0;	
	for (i = 0; i < xyzv.length/4; i++) {
		x = xyzv[4*i];
		y = xyzv[4*i+1];
		z = xyzv[4*i+2];
		w = xyzv[4*i+3];
		sw+= w;
		swx += w * x;
		swy += w * y;
		swz += w * z;
		swxx += w * x * x;
		swyy += w * y * y;
		swzz += w * z * z;
		swxy += w * x * y;
		swxz += w * x * z;
		swyz += w * y * z;
	}
	swx /= sw;
	swxx /= sw;
	swy /= sw;
	swyy /= sw;
	swz /= sw;
	swzz /= sw;
	swxy /= sw;
	swxz /= sw;
	swyz /= sw;
	swxx -= swx*swx;
	swyy -= swy*swy;
	swzz -= swz*swz;
	swxy -= swx*swy;
	swxz -= swx*swz;
	swyz -= swy*swz;
	return newArray(swxx,swxy,swxz,swyy,swyz,swzz);
}

function computeROIEntropy(xyzv) {
	/* Compute the entropy xyzv input with v wiewed as a distribution itself
	 *  Entropy the dentropy is defined by sum p log p when p is the intensity normlized by the 
	 *  sum of the intensity I / sum(I).
	 */
	norm = 0;
	for (i = 0; i < xyzv.length/4; i++) {
		maxi = norm += xyzv[4*i+3];		
	}
	entropy = 0;
	for (i = 0; i < xyzv.length/4; i++) {
		w = xyzv[4*i+3] / norm;
		if (w > 0.1) {
			entropy -= w * log(w);
		}
	}	
	return entropy;
}

function computeROITotalVolume(xyzv,threshold)  {	
	/* Compute the number of pixels value v which are above a threshold	 
	 *  Parameters
	 *  xyzv : 4xN array with coordinates and intensity
	 *  threshold : number
	 *  Result
	 *  volume (number)
	 */
	volume = 0;
	for (i = 0; i < xyzv.length/4; i++) {
		w = xyzv[4*i+3];
		if (w > threshold) {
			volume += 1.0;
		}
	}
	return volume;
}

function symeigen3x3(A) {
	/* Compute the eigen values of the symmetric 3x3 tensor A=[a11,a12,a13,a22,a23,a33]
	 *  https://arxiv.org/pdf/1306.6291.pdf
	 */
	A11 = A[0];
	A12 = A[1];
	A13 = A[2];
	A22 = A[3];
	A23 = A[4];
	A33 = A[5];
	b = A11+A22+A33;
	p = 0.5*((A11-A22)*(A11-A22)+(A11-A33)*(A11-A33)+(A22-A33)*(A22-A33)) + 3 * (A12*A12+A13*A13+A23*A23);	
	q = 18 * (A11*A22*A33+3*A12*A13*A23) + 2*(A11*A11*A11+A22*A22*A22+A33*A33*A33) + 9 * (A11+A22+A33)*(A12*A12+A13*A13+A23*A23) - 3 * (A11+A22)*(A11+A33)*(A22+A33) - 27*(A11*A23*A23+A22*A13*A13+A33*A12*A12);
	delta = acos(q/(2*sqrt(p*p*p)));
	print("d="+ q/(2*sqrt(p*p*p))) ;
	l1 = 1/3 * (b + 2*sqrt(p)*cos(delta/3));
	l2 = 1/3 * (b + 2*sqrt(p)*cos((delta+2*PI)/3));
	l3 = 1/3 * (b + 2*sqrt(p)*cos((delta-2*PI)/3));
	lambda = newArray(abs(l1),abs(l2),abs(l3));
	Array.sort(lambda);
	return lambda;	
}
