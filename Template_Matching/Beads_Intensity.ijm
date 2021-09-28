@Integer (label="Channel",description="Index of channel to process", value=1) channel
#@String (label="Channel names", description="Comma separated name of the channels", value="") channels_str
#@String (lebel="Mode", choices={"Correlation","Blobs","Skip"}) mode
#@Float (label="Downsampling", description="Downsampling factor forcomputing correlation", value=1) zoom 
#@Float (label="Threshold", description="values between -1 and 1", value=0.5) threshold
#@Float (label="Diameter [um]",description="Diameter of the beads", value=3) diameter
#@Float (label="Smoothness",desciption="Smoothness of the fitting", value=5) smoothness
#@Float (label="Bandwidth [px]", value=1) bandwidth

/*
 * Measure the intensity of large beads (~ 3 um) like objects 
 * 
 * Use template matching to locate beads and fit each location
 * using an active contour.
 * 
 * Usage
 * -----
 * - Open an image
 * - Add ROI to the ROI manager (press t after selecting the ROI)
 * - Open the macro in Fiji's script editor and press "Run"
 * - Adjust the parameters and press OK
 * 
 * If there is not opened image, a test image is generated.
 * 
 * If not enough similar structures are detected lower the treshold and/or the downsampling.
 * 
 * If a few ROI are incorrect, delete the problematic ROI, run again and untick the Find Similar option.
 * 
 * If the object are not fitted well enough, untick the Find Similar option, decrease the smoothness and run again
 * this will update the ROI.
 * 
 * 
 * Parameters
 * ----------
 * Channel : the channel to process
 * Channel names : comma separated list of the names of the channels 
 * Downsampling : the donsampling factor used for computing the correlation map
 * Threshold : the threshold [-1,1] to detect peaks in the correlation map
 * Smoothness : contour smoothness > 0,typically between 0.1 and 10
 * Bandwidth : width in pixels of the band in which the intensity is measured along the contours
 * 
 * Note
 * ----
 * The macro uses script parameters 
 *
 */

if (nImages==0) {
	print("\n >>>> Generating a test image <<<< \n");
	createTestImage();
	channel = 1;
	channels_str = "Red,Green";
	zoom = 2;
	threshold = 0.8;
	mode = "Correlation";
	smoothness = 5;
	diameter = 30;
	print("Reset zoom to "+zoom + " and threshold to "+threshold);
}

if (roiManager("count")==0 && matches(mode, "Correlation")) {
	exit("Please select some ROI and add them to the ROI manager or use another mode.");
}

starttime = getTime();
setBatchMode("hide");

print(" - Preprocessing channel " + channel + " at zoom " + zoom + "x");
id0 = getImageID();
run("Remove Overlay");
img = preprocess(channel);

// ROI detection
if (matches(mode,"Correlation")) {
	print(" - Finding candidate regions by normalized cross correlation");
	indices = correlateAllRoiFromManager(img, zoom, threshold);
}else if (matches(mode,"Blobs")) {
	print(" - Finding candidate regions by blob detection");
	detectBlobs(diameter);
}

// Fitting
// determine the typical area of the objects
print(" - Fitting " + roiManager("count") + " contours with smooth curves");
fitAllCurves(img, smoothness);

// Filtering
print(" - filtering "+roiManager("count")+" ROIs by area and circularity.");
a0 = PI/4 * diameter *diameter;
print("   >> Typical area :"+ a0+ " ["+0.75*a0+"-"+1.5*a0+"]");
filterROI(newArray("Circ.","Area"), newArray(0.9,0.75*a0), newArray(1,1.5*a0));

// Measurement
selectImage(id0);
print(" - Measuring mean intensity along "+roiManager("count")+" contours on image " + getTitle());
measure(id0, bandwidth, channels_str);
selectImage(img); close();
setBatchMode("exit and display");

// reset the display
run("Select None");
Stack.getDimensions(width, height, channels, slices, frames);
for (c=1;c<channels;c++) {
	Stack.setChannel(c);
	resetMinAndMax();
	run("Enhance Contrast", "saturated=0.2");
}
roiManager("show all without labels");
print("Elapsed time " + (getTime() - starttime)/1000 + " seconds.");

function createTestImage() {
	if (isOpen("ROI Manager")){selectWindow("ROI Manager");run("Close");}
	if (isOpen("Results")){selectWindow("Results");run("Close");}
	N = 500;
	newImage("Test Image", "32-bit composite-mode", N, N, 2, 1, 1);	
	D = 30;
	setColor(255);
	for (j = 0; j < 5; j++) {
		for (i = 0; i < 5; i++) {
			x = N*(i+1)/6.5 + D * random;
			y = N*(j+1)/6.5 + D * random;			
			run("Select None");
			makeOval(x, y, D, D);
			run("Draw", "stack");
			if (i==0 && j==0) {
				makeRectangle(x-5, y-5, D+10, D+10);
				roiManager("add");
				Roi.getCoordinates(xpoints, ypoints);
			}
		}
	}
	setColor(1);
	for (k = 0; k < 30; k++) {
		w = N/10*random;
		h = N/10*random;
		makeRectangle((N-w)*random, (N-h)*random, w, h);
		overlap = false;
		for (i = 0; i < xpoints.length && !overlap; i++) {
			if (Roi.contains(xpoints[i], ypoints[i])) {
				overlap = true;
			}
		}
		if (!overlap) {
			run("Draw", "stack");
		}
	}
	run("Select None");
	run("Gaussian Blur...", "sigma=2 stack");
	M = getValue("Max");
	run("Macro...", "code=v=v+"+(M/100)+"*(0.5+0.5*sin(0.1*x*d/w))");
	run("Add Specified Noise...", "standard="+M/20+" stack");
	Stack.setChannel(1); run("Enhance Contrast", "saturated=0.01");
	Stack.setChannel(2); run("Enhance Contrast", "saturated=0.01");
	run("8-bit");
}


function detectBlobs(D) {
	run("Select None");
	D = 50;
	id0 = getImageID();
	run("Duplicate...", "duplicate channels=1");
	run("32-bit");
	id1 = getImageID();
	run("Gaussian Blur...", "sigma="+D/4);
	run("Duplicate...", " ");
	id2 = getImageID();
	run("Gaussian Blur...", "sigma="+D);
	imageCalculator("Subtract 32-bit", id2, id1);
	selectImage(id1); close();
	selectImage(id2);
	getStatistics(area, mean, min, max, std, histogram);
	run("Find Maxima...", "prominence="+std/5+" output=[Point Selection]");
	Roi.getCoordinates(x, y);
	for (i=0;i<x.length;i++) {
		makeOval(x[i]-D/2, y[i]-D/2, D, D);
		roiManager("add");
	}
	close();
	selectImage(id0);
	roiManager("show all without labels");
}

function preprocess(channel) {
	run("Select None");
	run("Duplicate...", "title=preprocessed duplicate channels="+channel);
	run("32-bit");
	run("Gaussian Blur...", "sigma=2");
	run("Convolve...", "text1=[-1 -1 -1\n-1 8 -1\n-1 -1 -1\n] normalize");
	getStatistics(area, mean, min, max, std, histogram);
	run("Subtract...", "value="+std);
	run("Min...", "value=0");
	return getImageID();
}

function filterROI(key,lb,ub) {
	/*
	 * Filter ROI for which all keys values are above lb and lower than ub
	 */
	N = roiManager("count");
	for (i = N-1; i > 0; i--) {
		roiManager("select", i);
		ok = true;
		for (k = 0; k<key.length; k++) {
			v = getValue(key[k]);
			if (v < lb[k] || v > ub[k]) {
				ok = false;
				print(" - Remove ROI " + i + " with "+key[k]+"="+v+" <"+lb[k]+" or >"+ub[k]);
				roiManager("delete");	
				break;
			}
		}
	}
} 

function measure(id, bandwidth, channels_str) {
	selectImage(id);
	title = getTitle();
	n = roiManager("count");
	channel_names = split(channels_str,",");
	nr0 = nResults;
	for (i = 0; i < n; i++) {
		roiManager("select", i);
		setResult("Image", nr0 + i, title);
		setResult("Name", nr0 + i, Roi.getName());
		setResult("ID", nr0 + i, i+1);
		setResult("Area", nr0 + i, getValue("Area"));
		setResult("Circ.", nr0 + i, getValue("Circ."));
		vals = measureContourIntensity(id,bandwidth);
		for (c = 0; c < vals.length; c++) {
			if (c < channel_names.length) {
				setResult(channel_names[c], nr0 + i, vals[c]);
			} else {
				setResult("Channel "+(c+1), nr0 + i, vals[c]);
			}
		}
	}
	updateResults();
}

function measureContourIntensity(id,bandwidth) {
	// sum (img*mask*gaussian) / sum (mask*gaussian) for each channel
	tmax = 3 * bandwidth;
	run("Interpolate", "interval=0.5 smooth adjust");
	Roi.getSplineAnchors(x, y);
	dx = diff(x);
	dy = diff(y);
	L = x.length;
	K = Math.ceil(4*tmax);
	xt = newArray(K*L);
	yt = newArray(K*L);
	w = newArray(K*L);
	for (i=0;i<L;i++) {
		n = sqrt(dx[i]*dx[i]+dy[i]*dy[i]);
		for (k=0;k<K;k++) {
			t = -tmax + 2 * tmax * k / K;
			dt = t / bandwidth;
			xt[k+K*i] = x[i] + t * dy[i] / n;
			yt[k+K*i] = y[i] - t * dx[i] / n;
			w[k+K*i] = exp(-0.5*dt*dt);
		}
	}
	selectImage(id);
	Stack.getDimensions(width, height, channels, slices, frames);
	values = newArray(channels);
	for (c = 1; c <= channels; c++) {
		Stack.setPosition(c, 1, 1);
		s = 0;
		sw = 0;
		for (i = 0; i < K*L; i++) {
			s += getPixel(xt[i], yt[i]) * w[i];
			sw +=  w[i];
		}
		values[c-1] = s/sw;
	}
	return values;
}

////////////////////// Template Matching /////////////////////////////////////

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
		roiname = Roi.getName;
		getBoundingRect(template_x, template_y, template_width, template_height); 
		selectImage(small);
		makeRectangle(round(template_x/zoom), round(template_y/zoom), round(template_width/zoom), round(template_height/zoom)); 
		T = convertImageToArray(small);
		C = matchTemplateNormXCorr(I,T);
		xcorr = createImageFromArray("correlation","32-bit",C);
		run("Find Maxima...", "prominence="+threshold+" strict exclude output=[Point Selection]");
		Roi.getCoordinates(x, y);
		print("  >> Found " + x.length + " candidates for template #" + i + " '" + roiname + "'");
		for (j = 0; j < x.length; j++) {
			xj = zoom * x[j];
			yj = zoom * y[j];
			if (xj < image_width - template_width && yj < image_height - template_height) {
				dmin = 1e9;
				// add only roi away from existing one
				K = roiManager("count");
				for (k = 0; k < K; k++) {
					roiManager("select", k);
					getBoundingRect(xk, yk, width, height);
					dx = xk - xj;
					dy = yk - yj;
					dmin = minOf(dmin, dx*dx+dy*dy);
				}
				if (dmin > template_width * template_height) {
					selectImage(id);
					makeRectangle(xj, yj, template_width, template_height);
					roiManager("add");
					//print("  >> Adding " + Roi.getName + " distance to closest:" + sqrt(dmin));
				}
			}
		}
		selectImage(xcorr); close();
	}
	
	selectImage(small); close();
	
	N1 = roiManager("count");
	indices = newArray(N1-N+1);
	for (i = 0; i < indices.length; i++) {
		indices[i] = N + i;
	}
	return indices;
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


function createImageFromArray(title,type,A) {
	/* Create a 2D image from the array A
	 * The two first indices are the width and height
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



///////////////////// CURVE FITTING /////////////////////
function fitAllCurves(id,smoothness) {
	selectImage(id);
	n = roiManager("count");
	for (i = 0; i < n; i++) {
		showProgress(i, n);
		roiManager("select", i);
		run("Fit Circle");
		getSelectionBounds(x, y, width, height);
		fitCurve(500, width / 5, smoothness);
		roiManager("select", i);
		run("Restore Selection");
		roiManager("update");
	}
}

function fitCurve(niter,tmax,smooth) {
	getSelectionBounds(x, y, width, height);
	getPixelSize(unit, pw, ph);
	if (x!=0 && y!=0) {
		run("Interpolate", "interval="+5+"  smooth adjust");
		Roi.getSplineAnchors(x, y);
		delta = 1;
		for (i = 0; i < niter && delta > 0.01; i++) {
			delta = gradTheCurve(x,y,1+tmax*exp(-10*i/niter));
			fftSmooth(x,y,0.1+(smooth-0.1)*(1-exp(-10*i/niter)),true);
			Roi.setPolygonSplineAnchors(x, y);
			if (getValue("Area") < width * height * pw * ph / 5) {
				break;
			}
		}
		Roi.setPolygonSplineAnchors(x, y);
	}
}

function gradTheCurve(x,y,tmax) {
	// find the contour
	dx = diff(x);
	dy = diff(y);
	delta = 0; // displacement of the curve
	K = Math.ceil(3*tmax);
	for (i = 0; i < x.length; i++) {
		n = sqrt(dx[i]*dx[i]+dy[i]*dy[i]);
		swx = 0;
		sw = 0;
		for (k = 0; k < K; k++) {
			t = -tmax + 2 * tmax * k / (K-1);
			xt = x[i] + t * dy[i] / n;
			yt = y[i] - t * dx[i] / n;
			I = getPixel(xt,yt);
			w = I*exp(-0.5*(t*t)/(tmax*tmax));
			swx += w * t;
			sw += w;
		}
		if (sw>0) {
			tstar = swx / sw;
			vx = tstar * dy[i] / n;
			vy = - tstar * dx[i] / n;
			V = sqrt(vx*vx+vy*vy);
			x[i] = x[i] + vx / V;
			y[i] = y[i] + vy / V;
			delta = delta + V;
		}
	}
	return delta / x.length;
}

function diff(u) {
	// centered finite difference
	n = u.length;
	d = newArray(n);
	for (i = 1; i < n - 1; i++) {
		d[i] = u[i+1] - u[i-1];
	}
	d[n-1] = u[1] - u[n-1];
	d[0] = u[1] - u[n-1];
	return d; 
}


function fftSmooth(_real,_imag,_threshold,_hard) {
	//  Smoothing by Fourier coefficient hard thresholding
	N = nextPowerOfTwo(_real.length);
	real = Array.resample(_real,N);
	imag = Array.resample(_imag,N);
	fft(real,imag);
	maxi = 0;	
	for (i = 3; i < real.length; i++) {
		m = sqrt(real[i]*real[i]+imag[i]*imag[i])/_real.length;
		if (_hard) {
			if (m < _threshold) {
				real[i] = 0;
				imag[i] = 0;
			}
		} else {
			real[i] = (1-_threshold/m)*real[i];
			imag[i] = (1-_threshold/m)*imag[i];
		} 
		if (m>maxi) {maxi=m;}
	}
	ifft(real,imag);
	real = Array.resample(real,_real.length);
	imag = Array.resample(imag,_real.length);
	for (i = 0; i < real.length; i++) {
		_real[i] = real[i];
		_imag[i] = imag[i];
	}
}

function nextPowerOfTwo(_N) {
	// Return the next power of two integer
	return pow(2.0, 1+floor(log(_N)/log(2)));
}

function padarray(_src,_N) {
	//  Extend the size of the array to size N and fill it with zeros 
	if (_N > _src.length) {
	dest = newArray(_N);
	for (i = 0; i < _src.length; i++) {
		dest[i] = _src[i];
	}
	} else {
		dest=_src;
	}
	return dest;
	
}

function fft(_real,_imag) {
	// Forward Fast fourier transform
	_fft_helper(_real,_imag,0,_real.length,true);
}

function ifft(_real,_imag) {
	// Inverse fast fourier transform
	_fft_helper(_real,_imag,0,_real.length,false);
	for(k = 0; k < _real.length; k++) {
		_real[k] = _real[k] / _real.length;
		_imag[k] = _imag[k] / _imag.length;
	}
}

function _fft_helper(ax,ay,p,N,forward) {
	// Fast Fourier Transform helper function
    if(N < 2) {
    } else {
        separate(ax,p,N);
        separate(ay,p,N);        
        _fft_helper(ax, ay, p,       N/2, forward);  
        _fft_helper(ax, ay, p + N/2, N/2, forward);                               
        for(k = 0; k < N/2; k++) {
        	// even
            ex = ax[k    +p];
            ey = ay[k    +p];
            // odd
            ox = ax[k+N/2+p];
            oy = ay[k+N/2+p];
            // twiddle factor
            if (forward == true) {
            	wx = cos(-2*PI*k/N);
            	wy = sin(-2*PI*k/N);
            } else {
            	wx = cos(2*PI*k/N);
            	wy = sin(2*PI*k/N);
            }
            sx = wx * ox - wy * oy;
            sy = wx * oy + wy * ox;
            ax[k    +p] = ex + sx;
            ay[k    +p] = ey + sy;
            ax[k+N/2+p] = ex - sx;
            ay[k+N/2+p] = ey - sy;                      
        }           
    }
}

function separate(a, p, n) {
	// 	Separate odd and even
	b = newArray(n/2);	
    for(i=0; i < n/2; i++) {
        b[i] = a[2*i+1+p];        
    }
    for(i=0; i < n/2; i++) {
        a[i+p] = a[2*i+p];        
    }
    for(i=0; i < n/2; i++) {
        a[i+n/2+p] = b[i];        
    }
}