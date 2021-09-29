#@Integer (label="Channel",value=1) channel
#@Float (label="Downsampling", value=4) zoom 
#@Float (label="Threshold", value=0.5) threshold
#@Boolean (label="Average", value=true) average

/*
 * Template matching using normalized cross correlation
 * 
 * Select several ROIs in the image and press run to find similar objects using 
 * normalized cross correlation
 * 
 * 
 * Jerome Boulanger for Yue Liu
 */

if (nImages==0) {
	if (isOpen("ROI Manager")) {selectWindow("ROI Manager"); run("Close");}
	run("Dot Blot");
	makeRectangle(22, 28, 28, 28);	
	roiManager("add");	
}

setBatchMode("hide");
starttime = getTime();
id0 = getImageID();
run("Remove Overlay");
img = preprocess(channel);
indices = correlateAllRoiFromManager(img, zoom, threshold);
selectImage(img); close();
selectImage(id0);
setBatchMode("exit and display");
roiManager("show all without labels");

if (average) {
	setBatchMode(true);
	avg = averageAllROI(id0);
	setBatchMode(false);
}

print("Elapsed time " + (getTime() - starttime)/1000 + " seconds.");

function preprocess(channel) {
	run("Select None");
	run("Duplicate...", "duplicate channels="+channel);
	run("32-bit");
	run("Gaussian Blur...", "sigma=2");
	run("Convolve...", "text1=[-1 -1 -1\n-1 8 -1\n-1 -1 -1\n] normalize");
	return getImageID();
}

function correlateAllRoiFromManager(id,zoom,threshold) {
	/*
	 * Detect similar ROI by correlation
	 * Parameter
	 * ---------
	 * id : image identifier
	 * zoom : downsampling factor
	 * threshold : threshold of the normalized cross correlation [0,1]
	 * Result
	 * ------
	 * Detected ROI are added to the RI Manager
	 * Return the list of new ROI
	 */
	selectImage(id);
	run("Duplicate...", " ");
	image_width = getWidth();
	image_height = getHeight();
	run("Size...", "width="+round(image_width/zoom)+" height="+round(image_height/zoom)+" depth=1 constrain average interpolation=Bilinear");
	small = getImageID();
	I = convertImageToArray(small);
	N = roiManager("count");
	for (i = 0; i < N; i++) {
		selectImage(id);
		roiManager("select", i);
		getBoundingRect(template_x, template_y, template_width, template_height); 
		selectImage(small);
		makeRectangle(round(template_x/zoom), round(template_y/zoom), round(template_width/zoom), round(template_height/zoom)); 
		T = convertImageToArray(small);
		C = matchTemplateNormXCorr(I,T);
		xcorr = createImageFromArray("correlation","32-bit",C);
		run("Find Maxima...", "prominence="+threshold+" strict exclude output=[Point Selection]");
		Roi.getCoordinates(x, y);
		for (j = 0; j < x.length; j++) {
			if (zoom*x[j] >= image_width - template_width || zoom*y[j] >= image_height - template_height) {
				ok = false;
			} else {
				ok = true;
				// add only roi away from existing one
				K = roiManager("count");
				for (k = 0; k < K && ok; k++) {
					roiManager("select", k);
					getBoundingRect(xk, yk, width, height);
					dxk = xk - zoom * x[j];
					dyk = yk - zoom * y[j];
					if (dxk * dxk + dyk * dyk <= template_width * template_height) {
						ok = false;
					}
				}
				if (ok) {
					selectImage(id);
					makeRectangle(zoom*x[j], zoom*y[j], template_width, template_height);
					roiManager("add");
				}
			}
		}
		selectImage(xcorr); close();
	}
	selectImage(small); close();
	N0 = roiManager("count");
	for (i=N0-1; i>=0; i--) {
	}

	
	N1 = roiManager("count");
	indices = newArray(N1-N+1);
	for (i = 0; i < indices.length; i++) {
		indices[i] = N + i;
	}
	return indices;
}

function averageAllROI(id) {
	/*
	 * Average ROIs
	 * ROI are aligned at the top left corner
	 * The size of the average image is given by the largest ROI
	 * Parameter
	 * ---------
	 * id: image identifier
	 * Result
	 * ------
	 * returns the identifier of the average image
	 */
	N = roiManager("count");
	acc = 0;
	if (N > 0) {
		// find the size of the average image
		w = 0;
		h = 0;
		for (i = 0; i < N; i++) {
			roiManager("select", i);
			Roi.getBounds(x, y, width, height);
			w = maxOf(w, width); 
			h = maxOf(h, height);
		}
		newImage("Average", "32-bit black", w, h, 1);
		acc = getImageID();
		setPasteMode("Add");
		for (i = 0 ; i < N; i++) {
			selectImage(id);
			run("Select None");
			roiManager("select", i);
			run("Copy");
			selectImage(acc);
			run("Select None");
			run("Paste");
		}
		run("Select None");
		run("Divide...", "value="+N);
		run("Select None");
		run("Enhance Contrast", "saturated=0.01");
		}
	return acc;
}

function drawRectangles(x,y,w,h,zoom) {
	for (i=0;i<x.length;i++) {
		makeRectangle(x[i]*zoom, y[i]*zoom, w*zoom, h*zoom);
		Overlay.addSelection("yellow", 1);
	}
}

function matchTemplateNormXCorrXY(I,J,x,y) {
	/* Correlate the image I J with template J at position given by x,y
	 * Return the normalized cross correlation 
	 */
	wi = I[0];
	hi = I[1];
	wj = T[0];
	hj = T[1];
	C = newArray(x.length);
	for(i = 0; i<x.length;i++) {
		C[i] = nxcorr(I,J,x[i],y[i]);
	}
	return C;
}

function matchTemplateNormXCorr(I,J) {
	/* Correlate the image I J with template J
	 * Return the normalized cross correlation map
	 */
	wi = I[0];
	hi = I[1];
	C = newArray(2+wi*hi);
	C[0] = wi;
	C[1] = hi;
	for (yi = 0; yi < hi; yi++) {
		for (xi = 0; xi < wi; xi++) {
			C[2+xi+wi*yi] = nxcorr(I,J,xi,yi);;
		}
	}
	return C;
}

function nxcorr(I,J,xi,yi) {
	// Normalized cross correlation at xi,yi 
	wi = I[0];
	hi = I[1];
	wj = T[0];
	hj = T[1];
	sx = 0;
	sy = 0;
	sxx = 0;
	syy = 0;
	sxy = 0;
	n = 0;
	for (yj = 0; yj < hj; yj++) {
		for (xj = 0; xj < wj; xj++) {
			if (xi+xj < wi && yi+yj < hi) {
				vx = I[2+xi+xj+(wi)*(yi+yj)];
				vy = J[2+xj+(wj)*yj];
				sx += vx;
				sy += vy;
				sxx += vx*vx;
				syy += vy*vy;
				sxy += vx*vy;
				n++;
			}
		}
	}
	sx /= n;
	sy /= n;
	sxx /= n;
	syy /= n;
	sxy /= n;
	sxx -= sx*sx;
	syy -= sy*sy;
	sxy -= sx*sy;
	return sxy/sqrt(sxx*syy);
}

function matchTemplateSAD(I,J) {
	/* Correlate the image I J with template J
	 * Return the normalized cross correlation map
	 */
	wi = I[0];
	hi = I[1];
	wj = T[0];
	hj = T[1];
	C = newArray(2+wi*hi);
	C[0] = wi;
	C[1] = hi;
	for (yi = 0; yi < hi; yi++) {
		showProgress(yi,hi);
		for (xi = 0; xi < wi; xi++) {
			s = 0;
			n = 0;
			for (yj = 0; yj < hj; yj++) {
				for (xj = 0; xj < wj; xj++) {
					if (xi+xj < wi && yi+yj < hi) {
						vx = I[2+xi+xj+(wi)*(yi+yj)];
						vy = J[2+xj+(wj)*yj];
						s += abs(vx-vy);
						n++;
					}
				}
			}
			C[2+xi+wi*yi] = s / n;
		}
	}
	return C;
}


function computeImageMean(A){
	s = 0;
	for (i=2;i<A.length;i++) {s+=A[i];}
	return s / (A.length-2);
}


function createImageFromArray(title,type,A) {
	/* Create a 2D image from the array A
	 * The two first indices are the width and height
	 */ 
	width = A[0];
	height = A[1];
	if (A.length != 2+width*height) {
		print("Error in createImageFromArray: incompatible dimensions");
	}
	newImage(title, type, width, height, 1);
	for (y = 0; y < height; y++) {
		for (x = 0; x < width; x++) {
			setPixel(x,y,A[2 + x + width * y]);
		}
	}
	return getImageID();
}


function convertImageToArray(id) {
	/* Convert the current selection in the image to an array
	*/
	selectImage(id);
	getBoundingRect(x0, y0, width, height);
	A = newArray(2+width*height);
	A[0] = width;
	A[1] = height;
	for (y = 0; y < height; y++) {
		for (x = 0; x < width; x++) {
		v =  getPixel(x0+x, y0+y);
			A[2+x + width * y] = v;//getPixel(x0+x, y0+y);
		}
	}
	return A;
}

