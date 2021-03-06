#@Float(label = "aggregate threshold", value=2) aggregate_threshold
#@Float(label = "scan width", value=30) tmax
#@Float(label = "contour smoothness", value=1) smooth 
#@Integer(label = "maximum number of iterations", value=10) niter
#@Float(label = "bandwidth for measurement", value=10) bandwidth


/* Segment selected GUV and measure intensites
 * 
 * - Select ROIs on the first frame as close as possible to the GUV edge
 * - add the ROIs to the ROI manager
 * - Run the macro and set the parameter
 *
 * 	 scan width : how far (pixels) the contours is fitted 
 * 	 iteration : number of iterations
 * 	 bandwidth : band (pixels) for measurement (radius of around control points) 
 *   
 * The macro uses a contours based snake segmentation approach to find the
 * edge of the GUV based on the sum of normalized intensties of both channels. 
 * Raw intensities in each channels are then reported in a table for each frame.
 * 
 * 
 * Jerome Boulanger 2020-2021 for Saule Spokaite sspokaite@mrc-lmb.cam.ac.uk
 */


if (nImages==0) {
	print("Creating a test image");
	createTestImage(400,200,50);
	print("Number of ROI:" + roiManager("count"));
}

if (nSlices==1) {
	exit("The image is not a stack.");
}

if (isOpen("Results")) {selectWindow("Results");run("Close");}

// Check the number of ROI in the ROI manager
nroi = roiManager("count");
while (nroi == 0) {
	waitForUser("Please select ROI and add them to the ROI manager");
	nroi = roiManager("count");
}
print("Analysis of " + nroi + " regions.");

setBatchMode("hide");
trackAndMeasureIntensity(bandwidth, tmax, niter);
setBatchMode("exit and display");
print("Done");

function trackAndMeasureIntensity(bandwidth, tmax, niter) {
	// process all rois
	run("Select None");
	Overlay.remove;
	id0 = getImageID();
	print("Preprocessing image stack");
	id1  = combineChannels();
	selectImage(id0);
	mask = createMaskForAggregate(5);
	print("Analysing ROI");
	nroi = roiManager("count"); 
	for (roi = 0; roi < nroi; roi++) {
		processROI(id0,id1,mask,roi,tmax,niter,bandwidth);
	}
	selectImage(id1); close();
	selectImage(mask); close();
	Stack.getDimensions(width, height, nchannels, slices, frames);
	showGraphs(nchannels);
}

function processROI(id0,id1,mask,roi,tmax,niter,bandwidth) {
	// process a single roi
	selectImage(id0);
	Stack.getDimensions(width, height, channels, slices, frames);
	roiManager("select", roi);
	roiManager("add");
	current = roiManager("count") - 1;
	for (t = 1; t <= frames; t++) {
		trackGUV(id1,current,t,smooth);
		addOverlay(id0,mask,current,t,roi+1);
		measure(current,t,id0,mask,bandwidth,roi+1);
	}
	roiManager("select",current); roiManager("delete");
	selectImage(id0);
	Stack.setPosition(1,1,1);
}

function getPixelsValues(x,y) {
	// Return an array with the pixel intensity as position given by x and y
	a = newArray(x.length);
	for (i=0;i<x.length;i++) {
		a[i] = getPixel(x[i],y[i]);
	}
	return a;
}

function getPointsAndWeight(bandwidth) {
	tmax = 3 * bandwidth;
	run("Interpolate", "interval=0.5 smooth adjust");
	Roi.getSplineAnchors(x, y);
	dx = diff(x);
	dy = diff(y);
	K = Math.ceil(4*tmax);
	A = newArray(3+3*K*x.length);
	A[0] = 3;
	A[1] = K;
	A[2] = x.length;
	sw = 0;
	swx = 0;
	for (i=0;i<x.length;i++) {
		n = sqrt(dx[i]*dx[i]+dy[i]*dy[i]);
		for (k=0;k<K;k++){
			t = -tmax + 2 * tmax * k / K;
			xt = x[i] + t * dy[i] / n;
			yt = y[i] - t * dx[i] / n;
			dt = t / bandwidth;
			w = exp(-0.5*dt*dt);
			A[3+3*(k+K*i)] = xt;
			A[4+3*(k+K*i)] = yt;
			A[5+3*(k+K*i)] = w;
		}
	}
	return A;
}

function getPixelValuesInWeigtedBand(A) {
	K = A[1];
	L = A[2];
	I = newArray(K*L);
	for (i = 0; i < L; i++) {
		for (k = 0; k < K; k++){
			I[k+K*i] = getPixel(A[3+3*(k+K*i)],A[4+3*(k+K*i)]) * A[5+3*(k+K*i)];
		}
	}
	return I;
}


function getPixelValuesInBand(img,mask,bandwidth,frame) {
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
	m = newArray(K*L);
	selectImage(mask);
	Stack.setFrame(frame);
	for (i=0;i<L;i++) {
		n = sqrt(dx[i]*dx[i]+dy[i]*dy[i]);
		for (k=0;k<K;k++) {
			t = -tmax + 2 * tmax * k / K;
			dt = t / bandwidth;
			xt[k+K*i] = x[i] + t * dy[i] / n;
			yt[k+K*i] = y[i] - t * dx[i] / n;
			w[k+K*i] = exp(-0.5*dt*dt);
			m[k+K*i] = getPixel(xt[k+K*i],yt[k+K*i]);
		}
	}
	selectImage(img);
	Stack.getDimensions(width, height, channels, slices, frames);
	values = newArray(channels);
	for (c=1;c<=channels;c++) {
		Stack.setPosition(c, 1, frame);
		s = 0;
		sw = 0;
		for (i=0;i<K*L;i++) {
			s += getPixel(xt[i],yt[i]) * m[i] * w[i];
			sw +=  m[i] * w[i];
		}
		values[c-1] = s/sw;
	}
	return values;
}

function sum(x) {
	s = 0;
	for (i=0;i<x.length;i++) {s+=x[i];}
	return s;
}

function mean(x) {
	s = 0;
	for (i=0;i<x.length;i++) {s+=x[i];}
	return s/x.length;
}

function measure(roi,t,img,mask,bandwidth,roiindex) {
	// measure the intensity on a band defined by the line and the bandwidth
	// taking into account the mask. It also measures the intensity in a 10px band, 
	// 10px away from this line to define a background measure. The mask is used
	// here again.
	print("Measure roi:"+ roi + " at frame:" + t + " on  image " + getTitle());

	// background pixels
	selectImage(img);
	roiManager("select",roi);
	getPixelSize(unit, pixelWidth, pixelHeight);
	run("Enlarge...", "enlarge="+10*pixelWidth);
	run("Make Band...", "band="+10*pixelWidth);
	Roi.getContainedPoints(x0, y0);

	// measure mask
	selectImage(mask);
	roiManager("select",roi);
	Stack.setFrame(t);
	M0 = getPixelsValues(x0,y0);
	
	// measure in image
	values = getPixelValuesInBand(img,mask,bandwidth,t);
	
	selectImage(img);
	Stack.getDimensions(width, height, channels, slices, frames);
	setResult("Frame", t-1, t);
	for (c = 1; c <= channels; c++) {
		// intensity in the image
		roiManager("select",roi);
		Stack.setPosition(c, 1, t);
		setResult("ROI"+roiindex+"-ch"+c, t-1, values[c-1]);
		
		// background
		V0 = getPixelsValues(x0,y0);
		for (i=0;i<x0.length;i++) {V0[i] = V0[i]*M0[i];}
		setResult("BG"+roiindex+"-ch"+c, t-1, sum(V0)/sum(M0));
	}
}

function addOverlay(img,mask,roi,frame,label) {
	// add an overlay to the image
	print("Add overlay roi:"+ roi + " at frame:" + t);
	selectImage(mask);
	roiManager("select",roi);
	Roi.getCoordinates(x, y);	
	V = getPixelsValues(x,y);
	selectImage(img);
	Stack.setFrame(frame);
	for (i = 0; i < x.length-1; i++) {
		if ( V[i] > 0) {
			setColor("green");
		} else {
			setColor("red");
		}
		setLineWidth(2);
		Overlay.drawLine(x[i], y[i], x[i+1], y[i+1]);
		Overlay.add();
		Overlay.setPosition(1, 1, frame);
	}
	//Overlay.addSelection("red",1);
	Overlay.show();
	roiManager("select",roi);
	setColor("white");
	Overlay.drawString(""+label, mean(x), mean(y));
	Overlay.setPosition(1, 1, frame);
	Overlay.add();
	Overlay.show();
	
	// background
	selectImage(mask);
	getPixelSize(unit, pixelWidth, pixelHeight);
	run("Enlarge...", "enlarge="+10*pixelWidth);
	run("Make Band...", "band="+10*pixelWidth);
	roiManager("add"); roi1 = roiManager("count")-1;
	
	setThreshold(0.5, 1);
	run("Create Selection");
	if (getValue("Area") < getWidth*getHeight) {
		roiManager("add"); roi2 = roiManager("count")-1;

		roiManager("select", newArray(roi1,roi2));
		roiManager("and");
		roiManager("add"); roi3 = roiManager("count")-1;

		selectImage(img);
		roiManager("select", roi3);
		Overlay.addSelection("blue",1);
		roiManager("select", roi3);roiManager("delete");
		roiManager("select", roi2);roiManager("delete");
		roiManager("select", roi1);roiManager("delete");
	} else {
		selectImage(img);
		roiManager("select", roi1);
		Overlay.addSelection("blue",1);
		roiManager("select", roi1); roiManager("delete");
	}
	Overlay.setPosition(1, 1, frame);
	Overlay.add();

	Overlay.show();
}

function trackGUV(img,roi,frame,smooth) {
	print("Track  roi:"+ roi + " at frame:" + t);
	selectImage(img);
	roiManager("select", roi);
	getPixelSize(unit, pixelWidth, pixelHeight);
	run("Enlarge...", "enlarge="+(tmax/4*pixelWidth));
	run("Fit Circle");
	roiManager("Add");
	roiManager("select", roiManager("count")-2);
	roiManager("delete");
	roiManager("select", roiManager("count")-1);
	run("Interpolate", "interval="+5+"  smooth adjust");
	Roi.getSplineAnchors(x, y);
	Stack.setFrame(frame);
	delta = 1;
	for (i = 0; i < niter && delta > 0.5; i++) {
		delta = gradTheCurve(x,y,tmax);
		fftSmooth(x,y,smooth,true);
		Roi.setPolygonSplineAnchors(x, y);
	}
	Roi.setPolygonSplineAnchors(x, y);
	Roi.setPosition(1, 1, frame);
	roiManager("update");
}

function combineChannels() {
	// combine channels to extract the edge of the GUV
	run("Select None");
	id0 = getImageID();
	run("Duplicate...", "duplicate channels=1");
	id1 = getImageID();
	run("32-bit");
	run("Subtract Background...", "rolling=10 stack");
	Stack.getStatistics(voxelCount, avg, min, max, stdDev);
	run("Macro...", "code=v=(v-"+avg+")/"+stdDev+" stack");
	selectImage(id0);
	run("Duplicate...", "duplicate channels=2");
	id2 = getImageID();
	run("32-bit");
	Stack.getStatistics(voxelCount, avg, min, max, stdDev);
	run("Macro...", "code=v=(v-"+avg+")/"+stdDev+" stack");
	imageCalculator("Add stack",id1,id2);
	selectImage(id2); close();
	selectImage(id1);rename("Segmentation");
	run("Gaussian Blur...", "sigma=1 stack");
	run("Subtract Background...", "rolling=20 stack");
	run("Min...", "value=0 stack");
	resetMinAndMax();
	return id1;
}

function createMaskForAggregate(sz) {
	// Create a negative mask for large aggregates
	run("Duplicate...", "title=mask duplicate channels=1");
	id1 = getImageID();
	run("Median...", "radius="+sz+" stack");
	run("Minimum...", "radius="+sz+" stack");
	run("Maximum...", "radius="+sz+" stack");
	Stack.getStatistics(voxelCount, avg, min, max, stdDev);
	setThreshold(min, avg+aggregate_threshold*stdDev);
	run("Convert to Mask", "method=Default background=Default");
	run("Divide...", "value=255 stack");
	run("Enhance Contrast", "saturated=0.35");
	return id1;
}

function gradTheCurve(x,y,tmax) {
	// find the contour
	dx = diff(x);
	dy = diff(y);
	delta = 0; // displacement of the curve
	K = 4*tmax;
	for (i = 0; i < x.length; i++) {
		n = sqrt(dx[i]*dx[i]+dy[i]*dy[i]);
		swx = 0;
		sw = 0;
		for (k = 0; k < K; k++) {
			t = -tmax + 2 * tmax * k / K;
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
			x[i] = x[i] + vx;
			y[i] = y[i] + vy;
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

function createTestImage(W,H,T) {
	if(isOpen("ROI Manager")){selectWindow("ROI Manager");run("Close");}
	newImage("Test", "32-bit grayscale-mode", W, H, 2, 1, T);
	R  = 2;
	for (r = 0; r < R; r++) {
		x0 = (r+1) * W/(R+1);
		y0 = H / 2;
		d = H / 2;
		Stack.setPosition(1, 1, 1);
		makeOval(x0-(d+5)/2, y0-(d+5)/2,d+5,d+5);
		roiManager("add");
		for (t = 1; t <= T; t++) {
			x0 = x0 + random("gaussian");
			y0 = y0 + random("gaussian");
			for (c = 1; c <= 2; c++) { 
				Stack.setPosition(c, 1, t);
				if (c==1) {
					drawCircle(x0,y0,d,300);
				} else {
					drawCircle(x0,y0,d,(300+r*50)*(1/(1+exp(-10*(t-T/2)/T))));
				}
			}
		}
	}
	setBatchMode(true);
	run("Select None");
	newImage("Blotch", "32-bit grayscale-mode", W, H, 1, 1, T);
	run("Add Specified Noise...", "stack standard=5");
	makeRectangle(20, 20, W-40, H-40);
	run("Make Inverse");
	for (t=1;t<=T;t++){
		Stack.setPosition(c, 1, t);
		setColor(0);
		fill();
	}
	run("Select None");
	run("Gaussian Blur...", "sigma=20 stack");
	setThreshold(getValue("Mean")+3*getValue("StdDev"), 255);
	run("Convert to Mask", "method=Default background=Default");
	run("Gaussian Blur...", "sigma=5 stack");
	run("Duplicate...", "title=cpy duplicate");
	run("Merge Channels...", "c1=Blotch c2=cpy create");
	rename("Blotch");
	imageCalculator("Average stack", "Test", "Blotch");
	selectWindow("Blotch");close();
	run("Gaussian Blur...", "sigma=2 stack");
	run("Add Noise", "stack");
	for (c=1;c<=2;c++){
		Stack.setPosition(c, 1, T);
		resetMinAndMax();
	}
	run("8-bit");
	setBatchMode(false);
}

function drawCircle(x0,y0,d,val) {
	makeOval(x0-d/2, y0-d/2,d,d);
	run("Make Band...", "band=3");
	setColor(val);
	fill();	
}

function showGraphs(nchannels) {
	// Load the table and create graphs
	nroi = roiManager("count");
	Plot.create("Graph","Frame","Intensity")
	legend="";
	colors = getLutHexCodes("Ice",nchannels*nroi*2);
	mini=255;
	maxi=0;
	type = newArray("ROI","BG");
	for (i = 1; i <= nroi; i++) {
		for (c = 1; c <= nchannels; c++) {
			for (k = 0; k < type.length; k++) {
				T = Table.getColumn("Frame");
				I = Table.getColumn(""+type[k]+i+"-ch"+c);
				Array.getStatistics(I, min, max);
				mini = minOf(min,mini);
				maxi = maxOf(max,maxi);
				legend += ""+type[k]+i+"-ch"+c+"\t";
				Plot.setColor(colors[k+2*((c-1)+nchannels*(i-1))]);
				Plot.setLineWidth(2);
				Plot.add("line", T, I,""+type[k]+i+"-ch"+c);
			}
		}
	}
	Plot.setLimits(0, T[T.length-1], mini, maxi);
	Plot.setLegend(legend, "options");
}

function getLutHexCodes(name,n) {
	//setBatchMode("hide");
	newImage("stamp", "8-bit White", 10, 10, 1);
   	run(name);
   	getLut(r, g, b);
   	close();
   	r = Array.resample(r,n);
   	g = Array.resample(g,n);
   	b = Array.resample(b,n);
	dst = newArray(r.length);
	for (i = 0; i < r.length; i++) {
		dst[i] = "#"+dec2hex(r[i]) + dec2hex(g[i]) + dec2hex(b[i]);
	}   
	//setBatchMode("show");
	return dst;	
}

function dec2hex(n) {
	// convert a decimal number to its heximal representation
	codes = newArray("0","1","2","3","4","5","6","7","8","9","A","B","C","D","E","F");
	k1 = Math.floor(n/16);
	k2 = n - 16 * Math.floor(n/16);
	return codes[k1]+codes[k2]; 
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