#@File (label="Input File",style="open") inputfile
#@Integer (label="Number of frames [power of 2]", value=1024) N
#@Float (label="Distance to edge [um]", value=1.5) distance
#@Integer (label="Number of peaks", value=5) maxNFreq
#@Output Float fps
#@Output Float beat

/*
 * Cilia beating frequency analysis
 * 
 * Install stackreg (add update sites BIG-EPFL)
 * Run the macro as batch processing to process several files
 * 
 * The image sequence is stabilized using the pluging stackreg [1] twice with a rigid body
 * motion model: once on the full image and once on a subregion determined by segmenting 
 * the moving cell to account for residual rotations. The segmentation is performed by
 * thresholding the absolute value of the difference of the intensity with the mean 
 * of the image. A smoothing of the contour is then performed in fourier space to obtain 
 * a regular contour. A kymograph is then extracted from the contour and its 1d FFT along 
 * the temporal dimension computed [2]. Peaks in the max projection of the spectrum along the 
 * curvilinear abscissa enable to locate the main beat frequency. Kymograph, spectrum and graphs 
 * are stored for each processed file.
 * 
 * Parameters:
 * - Number of frames (N): limit the number of frames to analyse
 * - Distance to edge (distance): distance to the edge of the segmented cell where kymograph is computed
 * 
 * Output:
 * - Save the frequencies measures in a table "Beat Stats"
 * 
 * [1] P. Thévenaz, U.E. Ruttimann, M. Unser, "A Pyramid Approach to Subpixel Registration Based on Intensity," IEEE Transactions on Image Processing, vol. 7, no. 1, pp. 27-41, January 1998.
 * [2] C. M. Smith et al., ‘ciliaFA: a research tool for automated, high-throughput measurement of ciliary beat frequency using freely available software’, Cilia, vol. 1, no. 1, p. 14, Aug. 2012, doi: 10.1186/2046-2530-1-14.
 *
 * 
 * Jerome Boulanger 2020 for Girish Mali
 */
 
run("Close All");
setBatchMode("hide");
run("Set Measurements...", "  redirect=None decimal=3");

// Loading the image
print("Opening " + inputfile);
run("Bio-Formats Importer", "open=["+inputfile+"] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT use_virtual_stack");
name = getTitle;
id00 = getImageID();
dt = Stack.getFrameInterval();
run("Duplicate...", "duplicate range=1-"+N);
id0 = getImageID();

// Registration
print("Registration 1/2");
closeWindow("ROI Manager");
run("Select None");
run("StackReg", "transformation=[Rigid Body]");
Stack.setFrame(1);
if (createContourROI(0.5,distance)!=1){print("Error could not segment the image.");exit;};

print("Registration 2/2 (crop)");
run("Crop");
id1 = getImageID();
roiManager("Show None");
closeWindow("ROI Manager");
run("Select None");
run("StackReg", "transformation=[Rigid Body]");

// Segmentation
print("Segmentation");
Stack.setFrame(1);
if (createContourROI(0.5,distance)!=1){print("Error could not segment the image.");setBatchMode("exit and display");exit;};
roiManager("Save", replace(inputfile, ".nd2", "-rois.zip"));
roiManager("Show None");
roiManager("Select", 1);
run("Add Selection...");
run("Flatten","slice");
saveAsLog("Jpeg", replace(inputfile, ".nd2", "-img.jpg"));
close();

// Compute the kymograph
roiManager("select", 1);
run("Reslice [/]...", "output=1.000 start=Top avoid");
run("32-bit");
id2 = getImageID();

// preprocess the kymograph before 1D FFT
run("Duplicate...", "title=blur");
blur = getImageID();
run("Median...", "radius=5");
imageCalculator("Subtract 32-bit", id2, blur);
closeImage(blur);
selectImage(id2);
run("Enhance Contrast", "saturated=0.35");
saveAsLog("Tiff", replace(inputfile, ".nd2", "-kymo.tif"));

// Compute 1D fft
dft = compute1Dfft();
id3 = getImageID();
run("Enhance Contrast", "saturated=0.35");
saveAsLog("Tiff", replace(inputfile, ".nd2", "-dft.tif"));
close();

// post process 1d fft
dft = smooth1d(dft, 3);
dft = substractBaseline(dft, round(n/10));
k = newArray(dft.length);
for (i=0; i < k.length;i++) {
	k[i] = i / (2*dft.length*dt);
}
Plot.create("plot", "Frequency [Hz]", "log(|FFT|) [au]", k, dft);
Plot.update();
selectWindow("plot");
Plot.makeHighResolution("graph", 4);
selectWindow("graph");
saveAsLog("Jpeg", replace(inputfile, ".nd2", "-graph.jpg"));
close();
close();

// find peacks
freq = findPeakFreq(dft, dt);

fps=1/dt;
if (!isOpen("Beat Stats")) {
	Table.create("Beat Stats");
} else {
	selectWindow("Beat Stats");
}
n = Table.size();
Table.set("File", n, name);
Table.set("Frame rate [Hz]", n, fps);
beat = freq[0];
for (i = 0; i < minOf(maxNFreq, freq.length); i++) {
	colname = "Beat "+i+" [Hz]";
	Table.set(colname, n, freq[i]);
}
Table.update();

closeImage(id00);
closeImage(id0);
closeImage(id1);
closeImage(id2);
closeImage(id3);
closeWindow("ROI Manager");
setBatchMode("exit and display");
print("Done.");
print("");

function findPeakFreq(signal, dt) {
	n = 2 * signal.length; // the dft has only n/2 samples in ImageJ
	df = 1 / (n * dt);
	Array.getStatistics(signal, min, max, mean, stdDev);
	thres = 0.25*stdDev;
	print("Threshold fft: "+thres);
	maxima = Array.findMaxima(signal, thres);
	for (i = 0; i < maxima.length; i++) {
		maxima[i] = maxima[i] * df;
	}
	return maxima;
}

function substractBaseline(src,k) {
	n = src.length;
	dst = newArray(n);
	Array.getStatistics(src, min, max);
	for (i = 0; i < n; i++) {
		mini = max;
		for (j = i-k; j <= i+k; j++) {
			if (j>=0 && j < n) {
				if (src[j] < mini) {
					mini = src[j];	
				}
			}
		}
		dst[i] = src[i] - mini;
	}
	return dst;
}

function smooth1d(src,k) {
	n = src.length;
	dst = newArray(n);
	Array.getStatistics(src, min, max);
	for (i = 0; i < n; i++) {
		sum = 0;
		count = 0;
		for (j = i-k; j <= i+k; j++) {
			if (j >= 0 && j < n) {
				sum =  sum + src[i];
				count = count + 1;
			}
		}
		dst[i] = sum / count;
	}
	return dst;
}

function saveAsLog(fmt, name) {
	print("Saving "+ getTitle()+ " as "  + name);
	saveAs(fmt, name);
}

function compute1Dfft() {
	setBatchMode("hide");
	id = getImageID();
	nx = getWidth();
	ny = getHeight();
	for (x = 0; x <nx; x++) {
		selectImage(id);
		signal = newArray(ny);
		for (y = 0; y < ny; y++) {
			signal[y] = getPixel(x,y);
		}
		dft = Array.fourier(signal, "none");
		if (x==0) {
			nk = dft.length;
			acc = newArray(nk);
			newImage("DFT", "32-bit black", nx, nk, 1);
			dftid = getImageID();
		}
		selectImage(dftid);
		for (y = 0; y < nk; y++) {
			acc[y] = maxOf(acc[y], dft[y]);
			setPixel(x, y, log(dft[y])); 
		}
	}
	for (y = 0; y < nk; y++) {
		acc[y] = log(acc[y]);
	}
	selectImage(dftid); 
	run("Enhance Contrast", "saturated=0.35");
	setBatchMode("exit and display");
	return acc;
}

function createContourROI(s,distance) {
	clearROIManager();
	run("Duplicate...", "use");
	mask=getImageID();
	run("32-bit");
	mean = getValue("Mean");
	run("Subtract...", "value="+mean);
	run("Abs");
	run("Median...", "radius=3");
	threshold = s*getValue("StdDev");
	setThreshold(threshold, getValue("Max"));
	run("Convert to Mask");
	run("Median...", "radius=5");
	run("Fill Holes");
	run("Analyze Particles...", "size=1000-Infinity add");
	selectImage(mask);close();
	n = roiManager("count");
	if (n==1) {
		roiManager("select",0);
		run("Enlarge...", "enlarge="+distance);
		roiFftSmooth(0.25);
		run("Interpolate", "interval=2 smooth");
		run("Area to Line");
		roiManager("Add");
		roiManager("select", 1);
		run("To Bounding Box");
		run("Enlarge...", "enlarge=10");
		roiManager("Add");
	} else {
		print("error number of ROI = " + n);
	}
	return n;
}

function clearROIManager() {
	n = roiManager("count");
	if (n>0) {
		roiManager("select", Array.getSequence(n));
		roiManager("delete");
	}
}

function closeWindow(name){
	if (isOpen(name)){
		selectWindow(name);
		run("Close");
	}
}

function closeImage(id){
	if (isOpen(id)) {
		selectImage(id);
		close();
	}
}




function roiFftSmooth(_threshold) {
	run("Interpolate", "interval=1 smooth adjust");
	Roi.getCoordinates(x, y);		
	fftSmooth(x,y,_threshold);
	setColor("green");
	makeSelection("freehand", x, y);
}

/*
 * Smoothing by Fourier coefficient thresholding
 */
function fftSmooth(_real,_imag,_threshold) {
	N = nextPowerOfTwo(_real.length);
	real = Array.resample(_real,N);
	imag = Array.resample(_imag,N);
	fft(real,imag);	
	for (i = 2; i < real.length; i++) {
		m = sqrt(real[i]*real[i]+imag[i]*imag[i]);
		if (m < _threshold*_real.length) {
			real[i] = 0;
			imag[i] = 0;
		}
	}	
	ifft(real,imag);
	real = Array.resample(real,_real.length);
	imag = Array.resample(imag,_real.length);
	for (i = 0; i < real.length; i++) {
		_real[i] = real[i];
		_imag[i] = imag[i];
	}
}

/*
 * Return the next power of two integer
 */
function nextPowerOfTwo(_N) {
	return pow(2.0, 1+floor(log(_N)/log(2)));
}

/*
 * Extend the size of the array to size N and fill it with zeros 
 */
function padarray(_src,_N) {
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

/*
 * Forward Fast fourier transform
 */
function fft(_real,_imag) {
	_fft_helper(_real,_imag,0,_real.length,true);
}

/*
 * Inverse fast fourier transform
 */
function ifft(_real,_imag) {
	_fft_helper(_real,_imag,0,_real.length,false);
	for(k = 0; k < _real.length; k++) {
		_real[k] = _real[k] / _real.length;
		_imag[k] = _imag[k] / _imag.length;
	}
}

/*
 * Fast Fourier Transform helper function
 */
function _fft_helper(ax,ay,p,N,forward) {
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


/*
 * Separate odd and even
 */
function separate(a, p, n) {	
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

/*** END oF FFT **/