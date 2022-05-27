#@File(label="Input file") filename
#@String(label="Delimiter",value="-") del
#@Boolean(label="manual segmentation", value=false) manualsegmentation
#@String(label="Overlay", choices={"Curvature","DNE"},style="listBox") overlaytype
#@Boolean(label="Display info", value=true) iwantinfo
#@Boolean(label="Add colorbar", value=true) iwantcolorbar
#@Float(label="Colorbar min", value = -0.05) colorbarmin
#@Float(label="Colorbar max", value = 0.05) colorbarmax
#@Boolean(label="Save image as jpg and close", value=true) finalclear
#@Boolean(label="Save ROI", value=true) saveroi
#@Boolean(label="Used saved ROI", value=false) useSavedROI
#@Float(label="Marker scale", value=4) markerscale

/*
 * Quantify organoids morphology
 * 
 * usage: press run or batch
 * 
 * Set the parameters
 *  - set the filename
 *  - set the delimiter in the file name used to encode batch etc..
 *  - tick manual segmentation if you want to select manually the outline of the orgnanoid
 *  - set a threshold for segmentation (higher stricter)
 *  - set the minimum area to filter out smaller objects
 *  - tick the last option to save a jpg and close the image
 * 
 * Note:
 * - The input filename has the following format [Batch]-[Condition]-[Day]-[OrganoidIndex].extension
 * - If a roi file is present the roi will be use 
 * 
 * Output:
 * - Area/size (need to take into account pixel size)
 * - Shape measure (need to implement other methods)
 *   o fractal dimension
 *   o number of inflection points based on curvature
 *   o Dirichlet Normal Energy (DNE)
 * - Transparency measurement based on variations of intensity
 * 
 * https://onlinelibrary.wiley.com/doi/full/10.1002/ajpa.21489
 * https://arxiv.org/abs/1201.3097
 * 
 * Jerome Boulanger for Ilaria (2020)
 * 
 * #@File(label="Weka Segmentation model") wekamodel
 */

//
if (endsWith(filename,".tif")) {
	run("Close All");
	run("Set Measurements...", "  redirect=None decimal=10");
	closeWindow("ROI Manager");
	
	folder = File.getDirectory(filename);
	fname = File.getName(filename);
	tbl = "Organoid Morphology";
	n = initTable(tbl);
	processFile(n, folder, fname, tbl);
	wait(500);
	run("Collect Garbage");
}

function manualBatchProcess() {
	// process a folder containing tif files
	filelist = getFileList(folder);
	N = filelist.length;
	N = 1;
	for (i = 0; i < N; i++) {
		showProgress(i, N);
	    if (endsWith(filelist[i], ".tif")) { 
			processFile(i, folder, filelist[i], "Organoid Morphology");
	    } 
	}
}

function initTable(tblname) {
	// Initialize a result table if needed and return the next row index.
	if (!isOpen(tblname)) {
		Table.create(tblname);
	}
	selectWindow(tblname);
	return Table.size;
}

function processFile(idx, folder, fname, tblname) {
	// open and process file
	setBatchMode("hide");
	cond = parseConditionName(fname);
	print("Processing " + fname);
	roifilename = replace(folder+File.separator+fname,".tif",".roi");
	
	open(folder + File.separator + fname);
	id = getImageID();
	
	getPixelSize(unit, pixel_size, pixel_size);
	
	if (useSavedROI && File.exists(roifilename)) {
		roiManager("Open", roifilename);
	} else {	
		if (manualsegmentation) {
			setTool("freehand");
			setBatchMode("exit and display");
			waitForUser("Please draw the outline of the organoid");
			setBatchMode("hide");
			roiFftSmooth(0.1);
			roiManager("add");
			run("32-bit");
		} else {		
			run("32-bit");
			id = getImageID();
			print("segmenting image " + id);
			segment(id);
		}
	}
		
	//removeScaleBar();
	
	
	// Compute curvature
	roiManager("select", 0);	
	Roi.getSplineAnchors(x, y);
	run("Select None");
	
	
	length = getCurvilinearLength(x,y);
		
	// average rdius
	R0 = getAverageRadius(x,y,pixel_size);
	
	// curvature
	curvature = getContourCurvatureGeo(x,y,pixel_size);	
	E = log(getWeightedEnergy(curvature,length));

	// Inflection points
	inflex = getContourInflectionPoints(curvature, R0);
	K = inflex.length;

	// Dirichlet normal energy
	dnd = DND(x,y);
	DNE = log(getWeightedEnergy(dnd,length));
	
	// Fractal dimension
	D = roiFractalDimension(x,y,false);
	
	// Transparency
	T = transparency();

	I0 = measureIntensityOutside(id);
		
	//selectWindow(tblname);
	selectImage(id);
	roiManager("select", 0);
	Table.set("Batch", idx, parseBatchNumber(fname));
 	Table.set("Condition", idx, cond);
 	Table.set("Day", idx, parseDay(fname));
 	Table.set("Organoid", idx, parseOrganoidIndex(fname));
 	Table.set("Pixel size ["+unit+"]", idx, pixel_size);
 	
 	Table.set("Area", idx, getValue("Area"));
 	Table.set("Perimeter",idx,getValue("Perim."));
 	Table.set("Average radius", idx, R0);
 	Table.set("Roundness", idx, getValue("Round"));
 	Table.set("Aspect Ratio", idx, getValue("AR"));
 	Table.set("Feret", idx, getValue("Feret"));
 	Table.set("MinFeret", idx, getValue("MinFeret"));
 	Table.set("Minor", idx, getValue("Minor"));
 	Table.set("Major", idx, getValue("Major"));
 	Table.set("Circularity", idx, 100*getValue("Circ."));
 	Table.set("Inflections Points", idx, K);
 	Table.set("Average curvature", idx, E);
 	Table.set("DNE", idx, DNE);
 	Table.set("Fractal Dimension", idx, D);
 	Table.set("Mean Intensity", idx, getValue("Mean"));
 	Table.set("StdDev Intensity.", idx, getValue("StdDev"));
 	Table.set("Transparency", idx, T); 	
 	Array.getStatistics(dnd, min, max, mean, stdDev);
 	Table.set("Min DND", idx, min);
 	Table.set("Max DND", idx, max);
 	Array.getStatistics(curvature, Cmin, Cmax, Cmean, Cstd);
 	Table.set("Min curvature", idx, Cmin);
 	Table.set("Max curvature", idx, Cmax);
 	Table.set("Mean curvature", idx, Cmean);
 	Table.set("Std curvature", idx, Cstd); 	
 	Table.set("Intensity outside", idx, I0);
 	 	 
 	//roiManager("select", 0);
 	run("Select None");
 	// Overlays
	if (matches(overlaytype,"DNE")) {
		colorContour(x,y,dnd,markerscale);
		if (iwantcolorbar) {
			addColorBar(dnd);
		}
	} else {
		colorContour(x,y,curvature,markerscale);
		if (iwantcolorbar) {
			addColorBar(curvature);
		}
	}
	
	print("there are " + roiManager("count") + " rois");
	// drawing the inflection poinrs as red circles
	setColor("red");
	for (i = 0; i < K; i++) {
		Overlay.drawEllipse(x[inflex[i]]-markerscale, y[inflex[i]]-markerscale,2*markerscale, 2*markerscale);		
		Overlay.show;
		Overlay.setStrokeWidth(markerscale);
	}

	// Display information as text
	if (iwantinfo) {
		setColor("black");
		setFont("Fixed", 30, "bold");
		Overlay.drawString("TRANSPARENCY: "+round(T)+"\nCURVATURE: "+R0*Cstd+"\nINFLECTION POINTS: "+K+"\nRADIUS:"+R0, 0, 30);
		//Overlay.drawString("TRANSPARENCY: "+round(T)+"\nCURVATURE: "+E+"\nDNE: "+DNE+"\nINFLECTION POINTS: "+K+"\nFRACTAL DIM: "+D, 0, 30);
		Overlay.show;
	}
 	 	
 	selectImage(id);
 	run("Select None");
 	if (finalclear) { 
 		ofilename = replace(folder+File.separator+fname,".tif",".jpg");
 		print("Saving file \n" + ofilename);
 		saveAs("Jpeg", ofilename);
 		close(); 
 	}
 	print("there are " + roiManager("count") + " rois");
 	if (saveroi && roiManager("count") == 1) {  		
 		print("Saving file \n" + roifilename);
 		roiManager("select", 0);
 		roiManager("Save", roifilename);
 	} else {
 		print("There are more than one ROI.. cannot save it as .roi file.");
 	}
 	 	
 	run("Select None");
 	closeRoiManager();
 	run("Select None");
 	setBatchMode("exit and display");
}

function addColorBar(values) {
	// add a colorbar with min/max given by the array values.
	run("RGB Color");
	W = getWidth;
	H = getHeight;
	w = round(0.025*W);
	h = round(0.7*H);
	fs = 20; // font size
	img = getImageID();	
	newImage("Untitled", "8-bit white", w, h, 1);
	cb = getImageID();
	run("Macro...", "code=v=255*(h-y)/h");
	run("Ice");
	run("Select All");
	run("Copy");
	selectImage(img);
	makeRectangle(W-w-3*fs-10, (H-h)/2, w, h);
	run("Paste");
	run("Select All");
	if (colorbarmin == 0 && colorbarmax == 0) {
		Array.getStatistics(values, colorbarmin, colorbarmax, mean, stdDev);
	}
	nticks=5;
	for (i = 0; i < nticks; i++) {
		setFont("Arial", fs);
		setColor("black");
		str = d2s(colorbarmin+i/(nticks-1)*(colorbarmax-colorbarmin), 3);
		Overlay.drawString(str, W-3*fs-5, (H+h)/2-i/(nticks-1)*h);
	}
	Overlay.show();
	selectImage(cb);close();
	
}

function segment(id) {
	preprocess(id);	
	keepLargestROI();	
	fitContour(id);
	//roiManager("select", 0);
	//roiFftSmooth(0.1);
	roiManager("show none");
	run("Select None");
}

function preprocess(id) {	
	/* Preprocess the image (32-bit) and populate the roi manager */
			
	selectImage(id);
	w = getWidth();
	h = getHeight();
	
	// smooth the image
	run("Duplicate...","title=tmp");
	a = getImageID();
	run("Size...", "width="+(w/4)+" height="+(h/4)+" depth=1 constrain average interpolation=Bilinear");		
	for (i = 0; i < 10; i++) {
		run("Minimum...", "radius=2");
		run("Maximum...", "radius=2");
		run("Median...", "radius=2");
	}
	
	// create a background ans substract it
	selectImage(a);
	run("Duplicate...", "title=bg");
	b = getImageID();
	for (i=0;i<4;i++) {
		run("Maximum...", "radius=" + round(w/32));
		run("Gaussian Blur...", "sigma=10");
		run("Minimum...", "radius=" + round(w/32));
	}
	run("Gaussian Blur...", "sigma="+w/32);
	imageCalculator("Subtract", a, b);
	
	// rescale to the original size
	selectImage(a);
	run("Size...", "width="+w+" height="+h+" depth=1 constrain average interpolation=Bilinear");	
	
	// create ROI
	setAutoThreshold("Otsu");
	setOption("BlackBackground", false);
	run("Convert to Mask");
	run("Fill Holes");
	run("Analyze Particles...", "size=0-Infinity add");
	run("Select None");
	// clean up	images
	selectImage(a);close();
	selectImage(b);close();	
	selectImage(id);	
}

function fitContour(id) {
	/*Fit the contours of the image id, update the ROI manager*/	
	selectImage(id);
	// create a contour image	
	run("Select None");
	run("Duplicate...","title=tmp2");
	run("32-bit");
	run("Gaussian Blur...", "sigma=1");
	run("Median...", "radius=5");
	run("Gaussian Blur...", "sigma=1");
	run("Find Edges");
	id1 = getImageID();
	
	initROI(id);
	
	
	selectImage(id1);
	roiManager("select",0);
	run("Interpolate", "interval="+10+"  smooth adjust");
	Roi.getSplineAnchors(x, y);	
	delta=1;
	showStatus("Optimizing contour");
	niter = 200;
	for (iter = 0; iter < niter && delta > 0.1; iter++) {
		showProgress(iter, niter);
		delta = gradTheCurve(x,y,10+50*exp(-iter/10));	
		fftSmooth(x,y,1);				
	}
	roiManager("select",0);
	x = Array.reverse(x);
	y = Array.reverse(y);
	Roi.setPolygonSplineAnchors(x, y);		
	roiManager("update");
	selectImage(id1);close();
	selectImage(id);
	run("Select None");
	
}

function initROI(id) {	
	/* erode the roi and select the biggest resulting region*/
	selectImage(id);
	run("Select None");
	run("Duplicate...","title=tmp");	
	run("8-bit");	
	run("Select All");
	setColor(255);
	fill();	
	id2 = getImageID();
	roiManager("select",0);		
	run("Enlarge...", "enlarge=-80");
	run("Enlarge...", "enlarge=30");
	setColor(0);
	fill();
	setThreshold(0, 0);	
	run("Analyze Particles...", "size=10-Infinity add");
	roiManager("select", 0);
	roiManager("delete");
	keepLargestROI();
	roiManager("select", 0);	
	roiManager("update");
	run("Select None");
	selectImage(id2); close();	
	run("Select None");
}

function keepLargestROI() {
	/* keep the largest ROI from the ROI manager */
	n = roiManager("count");
	run("Select None");
	amax = 0;
	istar = -1;
	for (i = 0; i < n; i++) {
		roiManager("select",i);
		a = getValue("Area");
		if (a > amax) {
			istar = i;
			amax = a;
		}
	}
	for (i = n-1; i >= 0; i--) {		
		if (i != istar) {
			roiManager("select",i);
			roiManager("delete");
		}
	}	
	roiManager("select",0);
	
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
		if (sw > 0) {
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

function segmentWeka(model) {	
	/* Segment the image using a pixel classifier */
	run("Duplicate...", " ");
	run("Trainable Weka Segmentation");
	wait(1000);
	selectWindow("Trainable Weka Segmentation v3.3.2");
	call("trainableSegmentation.Weka_Segmentation.loadClassifier", model);
	call("trainableSegmentation.Weka_Segmentation.getResult");
	selectWindow("Classified image");
	selectWindow("Trainable Weka Segmentation v3.3.2");
	close();
	selectWindow("Classified image");
	setThreshold(level-0.1, level+0.1);
	run("Convert to Mask");
	
	run("Analyze Particles...", "size=0-Infinity add");
	// keep the largest ROI
	n = roiManager("count");
	amax = 0;
	istar = -1
	for (i = 0; i < n; i++) {
		roiManager("select",i);
		a = getValue("Area");
		if (a > amax) {
			istar = i;
			amax = a;
		}
	}
	for (i = n-1; i >= 0; i--) {
		roiManager("select",i);
		if (i != istar) {
			roiManager("delete");
		}
	}
	close();
	roiManager("select",0);
	roiFftSmooth(0.1);	
	selectImage(id0);
}

function segment_old() {
	// Segment the organoid
	run("Duplicate...","title=tmp");
	selectWindow("tmp");
	run("32-bit");
	// preprocessing on a downsampled image
	w = getWidth();
	h = getHeight();
	run("Size...", "width="+(w/8)+" height="+(h/8)+" depth=1 constrain average interpolation=Bilinear");
	for (i = 0; i < 5; i++) {
		run("Minimum...", "radius=1");
		run("Maximum...", "radius=1");
		run("Median...", "radius=1");
	}
	run("Subtract Background...", "rolling=100 light sliding");
	run("Size...", "width="+w+" height="+h+" depth=1 constrain average interpolation=Bilinear");
	// Otsu thresholding
	getStatistics(area, mean, min, max, std, histogram);
	setThreshold(min, mean-alpha*std);
	run("Convert to Mask");
	run("Fill Holes");
	//makeOval(getWidth/8, getHeight()/8, 6/8*getWidth(), 6/8*getHeight());
	run("Analyze Particles...", "size=0-Infinity add");
	roiManager("select",0);
	roiFftSmooth(0.1);
	roiManager("update");
	close();
}

function substractBackground() {
	w = getWidth();
	h = getHeight();
	run("Size...", "width="+(w/8)+" height="+(h/8)+" depth=1 constrain average interpolation=Bilinear");
	for (i = 0; i < 5; i++) {
		run("Minimum...", "radius=1");
		run("Maximum...", "radius=1");
		run("Median...", "radius=1");
	}
	run("Subtract Background...", "rolling=100 light sliding");
	run("Size...", "width="+w+" height="+h+" depth=1 constrain average interpolation=Bilinear");
}

function measureIntensityOutside(id) {
	/* Measure intensity outside the organoid
	Use a simple threshold to discard the darker organoid and shrink the selection
	to avoid taking into account regions from inside the organoid. return the mean intensity
	*/
	selectImage(id);	
	setAutoThreshold("Minimum dark");
	run("Create Selection");
	run("Enlarge...", "enlarge=-50");
	val = getValue("Mean");
	run("Select None");	
	resetThreshold;
	return val;
}

function transparency() {
	// Measure the transparency as the Mean of a filtered image
	run("Duplicate...", " ");
	run("32-bit");
	run("Gaussian Blur...", "sigma=2");
	run("Convolve...", "text1=[-1 -1 -1\n-1 8.2 -1\n-1 -1 -1\n] normalize");
	v = getValue("Mean");
	close();
	return v;
}

function parseBatchNumber(fname) {
	dst = split(fname, del);
	return parseInt(substring(dst[0], 1));
}

function parseConditionName(fname) {
	dst = split(fname, del);
	return dst[1];
}

function parseDay(fname) {
	dst = split(fname, del);
	return parseInt(substring(dst[2], 3));
}

function parseOrganoidIndex(fname) {
	dst = split(fname, del);
	dst = split(dst[3],'.');
	return parseInt(dst[0]);
}

function removeScaleBar(){
	a0 = getValue("Area");
	setThreshold(255, 255);	
	run("Create Selection");
	if (getValue("Area")!=a0) {
		if (getValue("Area") != getWidth() * getHeight()) {
			run("Enlarge...", "enlarge=2");
			run("Median...", "radius=10");
			run("Select None");
		}
	}
}

function getContourInflectionPoints(curvature,r0) {
	// return list of index where the curvature is more than 2 time the minor axis.
	//r0 = getValue("Major");
	//tmp = smooth1d(curvature,2);
	tmp = curvature;
	//Array.getStatistics(tmp, min, max, mean, stdDev);	
	threshold = 0.5/r0;//mean + 0.5*stdDev;		
	inflx = findPeaks(tmp, threshold, 10);	
	
	/*
	xp = newArray(inflx.length);
	yp = newArray(inflx.length);
	for (i=0;i<inflx.length;i++){
		xp[i] = inflx[i];
		yp[i] = tmp[inflx[i]]; 
	}
	
	id = getImageID();
	Plot.create("Title", "X","Y");
	Plot.add("line", Array.getSequence(curvature.length), curvature);
	Plot.setColor("green");
	Plot.add("line", Array.getSequence(curvature.length), tmp);	
	xl = newArray(0,curvature.length-1);
	yl = newArray(threshold,threshold);
	Plot.setColor("red");
	Plot.add("line",xl,yl);
	Plot.add("circle",xp,yp);
	Plot.update();
	selectImage(id);
	*/
	return inflx;
}

function findPeaks(A,t,p) {
	/* find local maxima in neighborhood [-p,p] above threshold */
	maximas = newArray(A.length);
	k = 0;
	for (i = 0; i < A.length; i++) if (A[i] > t) {			
		flag = true;
		for (j = i-p; j <= i+p && flag; j++) {			
			if (A[getCircBB(A.length,j)] > A[i]) {
				flag = false;				
			}
		}
		if (flag) {
			maximas[k] = i;
			k++;
		}
	}	
	return Array.trim(maximas, k);
}

function getCircBB(n,i) {
	// index with circular boundary conditions
	if (i >= 0 && i < n) return i;
	if (i < 0) return n+i;
	if (i >=  n) return i-n;
}

function getCurvilinearLength(x,y) {
	// return the distance between consecutive points
	s = newArray(x.length);
	for (i = 0; i < n-1; i++) {
		dx = x[i+1] - x[i];
		dy = y[i+1] - y[i];
		s[i] = sqrt(dx*dx+dy*dy);
	}
	dx = x[0] - x[x.length-1];
	dy = y[0] - y[y.length-1];
	s[x.length-1] = sqrt(dx*dx+dy*dy);
	return s;
}

function getWeightedEnergy(value,weight) {
	// compute the energy of the curvature (sum weigth * value^2) / sum(weight)
	e = 0;
	p = 0;
	for (i = 0; i < value.length; i++) {
		e = e + value[i] * value[i] * weight[i];
		p = p + weight[i];
	}
	return e / p;
}

function DNDhelper(x0,y0,x1,y1,x2,y2) {
	// Norm of the variation of the normal projected on the tangent times 
	// compute derivatives at point p-1 and p+1
	dx1 = x1 - x0;
	dx2 = x2 - x1;
	dy1 = y1 - y0;
	dy2 = y2 - y1;
	s1 = sqrt(dx1 * dx1 + dy1 * dy1);
	s2 = sqrt(dx2 * dx2 + dy2 * dy2);
	dx1 = dx1 / s1;
	dx2 = dx2 / s2;
	dy1 = dy1 / s1;
	dy2 = dy2 / s2;
	// derivative of the normal n=(dy,-dx) 
	dnx = (dy2-dy1)/s1;
	dny = -(dx2-dx1)/s1;
	// project on tangent u=(dx,dy): dn.u
	dnu = dnx * dx1 + dny * dy2;
	return sqrt(dnu*dnu);
}

function DND(x,y) {
	n = x.length;
	e = newArray(n);
	for (i = 1; i < n-1; i++) {
		e[i] = DNDhelper(x[i-1],y[i-1],x[i],y[i],x[i+1],y[i+1]);
	}
	e[0] = DNDhelper(x[n-1],y[n-1],x[0],y[0],x[1],y[1]);
	e[n-1] = DNDhelper(x[n-2],y[n-2],x[n-1],y[n-1],x[0],y[0]);
	return e;
}

function colorContour(x,y,value,width) {
	// color the contour with value as overlay
	lut = getLutHexCodes("Ice");	
	//Roi.getCoordinates(x, y);
	//Roi.getSplineAnchors(x, y);	
	if (colorbarmin==0 && colorbarmax==0) {
		Array.getStatistics(value, colorbarmin, colorbarmax, mean, stdDev);
		colorbarmin = - maxOf(abs(colorbarmax), abs(colorbarmin));
		colorbarmax = -colorbarmin;
	}
	for (i = 0; i < x.length-1; i++) {
		cidx = Math.floor((lut.length-1)*(value[i] - colorbarmin) / (colorbarmax-colorbarmin));
		cidx = minOf(lut.length-1, maxOf(0, cidx));
		setColor(lut[cidx]);
		Overlay.drawLine(x[i], y[i], x[i+1], y[i+1]);
		Overlay.show;
		Overlay.setStrokeWidth(width);
		//Overlay.setStrokeColor(lut[cidx]);
		//Overlay.addSelection(lut[cidx], 2);
	}
	// close the  contour
	cidx = Math.floor((lut.length-1)*(value[i] - colorbarmin) / (colorbarmax-colorbarmin));
	setColor(lut[cidx]);
	Overlay.drawLine(x[x.length-1], y[y.length-1], x[0], y[0]);
	Overlay.show;
	run("Labels...", "color=white font=12");
}

function getLutHexCodes(name) {
	//setBatchMode("hide");
	newImage("stamp", "8-bit White", 10, 10, 1);
   	run(name);
   	getLut(r, g, b);
   	close();
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


function getAverageRadius(x,y,pixel_size) {
	/* compute the average radius of the polygone defined by x,y */
	n = x.length;
	xbar = 0;
	ybar = 0;
	for (i = 0; i < n; i++) {
		xbar += x[i];
		ybar += y[i];
	}
	xbar /= n;
	ybar /= n;
	r = 0;
	for (i = 0; i < n; i++) {
		dx = x[i]-xbar;
		dy = y[i]-ybar;
		r += sqrt(dx*dx+dy*dy);
	}		
	return r/n * pixel_size;
}

function getContourCurvatureGeo(x,y,pixel_size) {
	// Compute cuvature with periodic boundary conditionsusing a geometric method
	n = x.length;
	kappa = newArray(n);
	for (i = 0; i < n; i++) {
		ip = getCircBB(n,i-1);
		in = getCircBB(n,i+1);		
		kappa[i] = curvHelper(x[ip],y[ip],x[i],y[i],x[in],y[in]) / pixel_size;		
	}
	return kappa;
}

function curvHelper(x0,y0,x1,y1,x2,y2) {
	// Use Heron formula to compute curvature
	// as the radius of the osculating circle
	// using a geometric approach
	//  https://scholar.rose-hulman.edu/cgi/viewcontent.cgi?article=1233&context=rhumj
	xa = x0 - x1;
	ya = y0 - y1;
	
	xb = x2 - x1;
	yb = y2 - y1;
	
	xc = x2 - x0;
	yc = y2 - y0;
	
	a = sqrt(xa*xa+ya*ya);
	b = sqrt(xb*xb+yb*yb);
	c = sqrt(xc*xc+yc*yc);
	//S = sqrt((a+b+c)*(-a+b+c)*(a-b+c)*(a+b-c)) / 4;
	//R = a*b*c/(2*S);	
	Delta = (xb*ya-xa*yb) / 2.0;
	Kappa = 4.0 * Delta / (a*b*c);
	return Kappa;	
}


function getContourCurvature(x,y) {
	// Compute cuvature with periodic boundary conditions	
	//x = smooth1d(x,2);
	//y = smooth1d(y,2);
	n = x.length;
	gamma = newArray(n);
	for (i = 1; i < n-1; i++) {
		dx = (x[i+1] - x[i-1])/2;
		dy = (y[i+1] - y[i-1])/2;
		dxx = x[i-1] - 2 * x[i] + x[i+1];
		dyy = y[i-1] - 2 * y[i] + y[i+1];
		gamma[i] = (dx * dyy - dy * dxx) / pow(dx*dx + dy*dy, 3.0/2.0);
	}
	// at 0 with perdiodic bc
	dx = x[1] - x[n-1];
	dy = y[1] - y[n-1];
	dxx = - x[n-1] + 2 * x[0] - x[1];
	dyy = - y[n-1] + 2 * y[0] - y[1];
	gamma[0] = (dx * dyy - dy * dxx) / pow(dx*dx + dy*dy, 3.0/2.0);
	// at n-1 with periodic bc
	dx = x[n-2] - x[0];
	dy = y[n-2] - y[0];
	dxx = - x[n-2] + 2 * x[n-1] - x[0];
	dyy = - y[n-2] + 2 * y[n-1] - y[0];
	gamma[n-1] = (dx * dyy - dy * dxx) / pow(dx*dx + dy*dy, 3.0/2.0);
	return gamma;
}


function diff(u) {
	n = u.length;
	d = newArray(n);
	for (i = 0; i < n - 1; i++) {
		d[i] = u[i+1] - u[i];
	}
	d[n-1] = u[0] - u[n-1];
	return d; 
}

function diff2(u) {
	n = u.length;
	d = newArray(n);
	for (i = 1; i < n-1; i++) {
		d[i] = - u[i-1] + 2 * u[i] - u[i+1];
	}
	d[n-1] = -u[n-2] + 2 * u[n-1] - u[0];
	d[0]   = -u[n-1] + 2 * u[0]   - u[1];  
	return d; 
}

function roiFractalDimension(x,y,doplot) {
	// Compute the fractal dimension of the contour
	// should be between 1 and 2..
	// https://arxiv.org/abs/1201.3097
	//run("Interpolate", "interval=1");
	//Roi.getCoordinates(x, y);	
	if (x.length > 10) {	
		if (0) {
		Array.getStatistics(x, xmin, xmax, xmean, xstdDev);
		Array.getStatistics(y, ymin, ymax, ymean, ystdDev);	
		for (i = 0; i < x.length; i++) {
			x[i]=(x[i]-xmean)/xstdDev;
			y[i]=(y[i]-ymean)/ystdDev;
		}
	}
	if (doplot==true) {
		Plot.create("curve","x","y",x,y);
	}
	m = powerspectrum(x, y);
	
	f = newArray(m.length-1);
	mf = newArray(m.length-1);
	for (i = 0; i <f.length; i++) {	
		f[i] = log(i+1);
		mf[i] = log(m[i+1]);
	}
	if (doplot==true) {
		Plot.create("curve","x","y",f,mf);
	}
	Fit.doFit("Straight Line", f, mf);
	if (doplot==true) {
		Fit.plot;
	}
	beta = Fit.p(1);
	return 1-(beta + 2) / 6;
	}
	return 1.0;
}

/*
 * Log powerspectrum
 */
function powerspectrum(_real,_imag) {
	//print(_real.length);	
	N = nextPowerOfTwo(_real.length);
	//print(N);
	real = Array.resample(_real,N);
	imag = Array.resample(_imag,N);
	//real = padarray(_real,N);
	//imag = padarray(_imag,N);
	fft(real,imag);
	m = newArray(real.length/2);
	for (i = 0; i < m.length; i++) {
		m[i] = real[i]*real[i]+imag[i]*imag[i];
	}
	return m;
}

function roiFftSmooth(_threshold) {
	// Smooth the contour using thresholding of the FFT
	run("Interpolate", "interval=1 smooth adjust");
	Roi.getCoordinates(x, y);		
	fftSmooth(x,y,_threshold);
	setColor("green");
	makeSelection("freehand", x, y);
	run("Interpolate", "interval=1 smooth adjust");
}

function fftSmooth(_real,_imag,_threshold) {
	//  Smoothing by Fourier coefficient hard thresholding
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
	//Plot.create("","","",_real,_imag);
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

/*** END oF FFT **/
function smooth1d(src,sigma) {
	n = src.length;
	dst = newArray(n);	
	k = Math.ceil(3*sigma);
	for (i = 0; i < n; i++) {
		sum = 0;
		count = 0;
		for (j = i-k; j <= i+k; j++) {
			d = (i-j)/sigma;
			w = exp(-0.5*d*d);
			sum =  sum + w * src[getCircBB(n,j)];
			count = count + w;			
		}
		dst[i] = sum / count;
	}
	return dst;
}

function closeRoiManager() {
	while (roiManager("count")>0) {roiManager("select", roiManager("count")-1);roiManager("delete");}
	closeWindow("ROI Manager");
	run("Select None");
}


function closeWindow(name){
	if (isOpen(name)){
		selectWindow(name);
		run("Close");
	}
}