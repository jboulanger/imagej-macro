#@String(label="mode",choices={"Track","Motion","Interpolate ROI", "Correct"}) mode
/*
 * Estimate translation using a bright reference point that is tracked over time or motion estimation
 * 
 * - Selecting a (rectangular) region of the image and run with the Track or Motion mode or add a set of 
 *   RI to the ROI manager to interpolate between them. ROI can be unsorted.
 * - The macro will save displements in a table
 * - Run again with the Correct mode to apply the estimated displacements
 * 
 * Jerome Boulanger 2021 for Florence update for Kashish
 */

if (nImages==0) {
	createTestImage();
	print("Run again the macro with track mode.");
	exit();
}

setBatchMode(true);
name = getTitle;
if (matches(mode, "Motion")) {	
	estimateMotion();
	print("Run again the macro with 'correct' mode.");
} else if (matches(mode, "Track")) {
	trackBrightSpot();
	print("Run again the macro with 'correct' mode.");
} else if (matches(mode, "Interpolate ROI")) {
	interpolateROI();
	print("Run again the macro with 'correct' mode.");
} else {
	run("Select None");
	run("Duplicate...", "title=[Stabilized+"+name+"] duplicate");
	correct();
}
print("Done");
setBatchMode(false);

function correct() {
	// Correct the drift fromt he stack using the value in the table
	Stack.getDimensions(width, height, channels, slices, frames);
	for (k = 2; k <= frames; k++){
		showProgress(k,frames);
		vx = getResult("X", k-1);
		vy = getResult("Y", k-1);			
		if (vx*vx+vy*vy>0.001) {
			Stack.setFrame(k);
			run("Translate...", "x="+-vx+" y="+-vy+" interpolation=Bicubic slice");
		}
	}
}

function interpolateROI() {
	/*
	 * Interpolate ROI poition from the ROI manager
	 * ROI do not need to be in order
	 */
	print("\\Clear");	
	Stack.getDimensions(width, height, channels, slices, frames);
	print("frames "+ frames);
	// compute the first and last frame of the ROI
	n = roiManager("count");
	T = newArray(n);
	X = newArray(n);
	Y = newArray(n);
	getPixelSize(unit, dx, dy);
	for (i = 0; i < n; i++) {
		roiManager("select", i);		
		Stack.getPosition(channel, slice, frame);
		T[i] = frame;
		if (matches(Roi.getType,"point")) {
			Roi.getCoordinates(xpoints, ypoints);
			X[i] = xpoints[0];
			Y[i] = ypoints[0];
		} else {			
			X[i] = getValue("X") / dx;
			Y[i] = getValue("Y") / dy;
		}		
	}
	Array.sort(T,X,Y);
	t0 = T[0];
	t1 = T[n-1];	
	for (t = 1; t <= frames; t++) {
		r = t - 1;		
		if (t <= t0) {			
			setResult("X", r, 0);
			setResult("Y", r, 0);
		} else if (t >= t1) {			
			setResult("X", r, X[n-1] - X[0]);
			setResult("Y", r, Y[n-1] - Y[0]);
		} else {	
					
			for (i = 0; i < n-2 && T[i] <= t; i++) {}	
			
			xa = X[i-1];
			ya = Y[i-1];
			ta = T[i-1];
			
			xb = X[i];
			yb = Y[i];			
			tb = T[i];
			
			print("t:"+t+":",i,", ta:"+ta+", tb:"+tb);
			setResult("X", r, xa + (t - ta) * (xb - xa) / (tb - ta) - X[0]);
			setResult("Y", r, ya + (t - ta) * (yb - ya) / (tb - ta) - Y[0]);
		}
		
	}
	updateResults();
}


function trackBrightSpot() {
	// estimate the drift by tracking a bright over time	
	run("Duplicate...", "title=tmp duplicate");
	run("Gaussian Blur 3D...", "x=0.5 y=0.5 z=0.1");	
	run("Subtract Background...", "rolling=10 stack");	
	
	Stack.getDimensions(width, height, channels, slices, frames);
	getPixelSize(unit, pixelWidth, pixelHeight);
	x = newArray(frames);
	y = newArray(frames);
	for (frame = 1; frame <= frames; frame++) {
		showProgress(frame,frames);
		Stack.setFrame(frame);
		run("Select None");
		getStatistics(area, mean, min, max, std, histogram);
		setThreshold(mean+2*std,max);
		run("Analyze Particles...", "size=1-Infinity pixel add slice");
		if (roiManager("count")>0) {
			istar = 0;
			if (frame==1) {
				a = 0;
				for (i=0;i<roiManager("count");i++) {
					roiManager("select", i);
					if (getValue("RawIntDen")>a) {
						a = getValue("RawIntDen");
						istar = i;
					}
				}
			} else {
				a = 1e9;
				for (i = 0; i < roiManager("count"); i++) {
					roiManager("select", i);
					dx = getValue("XM") / pixelWidth - x[frame-2];
					dy = getValue("YM") / pixelWidth - y[frame-2];
					d = dx*dx+dy*dy;					
					if (d < a) {
						a = d;
						istar = i;
					}
				}
			}			
			roiManager("select", istar);
			x[frame-1] = getValue("XM") / pixelWidth;
			y[frame-1] = getValue("YM") / pixelHeight;			
			Overlay.drawEllipse(x[frame-1]-2, y[frame-1]-2, 4,4);			
			Overlay.setPosition(frame);
			Overlay.add();
			Overlay.show();
			for (i = roiManager("count")-1; i>=0; i--) {
				roiManager("select", i);
				roiManager("delete");
			}			
		}
	}
	setResult("Frame",0,1);
	setResult("X",0,0);
	setResult("Y",0,0);
	for (i=1;i<x.length;i++) {
		setResult("Frame",i,i+1);
		setResult("X",i,x[i]-x[0]);
		setResult("Y",i,y[i]-y[0]);				
	}	
	close();
}


function estimateMotion() {
	// Estimate motion on the selected region
	print("\\Clear");
	run("Duplicate...", "title=[Stabilized] duplicate");	
	run("32-bit");
	run("Gaussian Blur 3D...", "x=0.5 y=0.5 z=1");	
	run("Subtract Background...", "rolling=10 stack");		
	
	id = getImageID();
	Stack.getDimensions(width, height, channels, slices, frames);
	
	for (k = 0; k < frames; k++){setResult("Frame", k, k+1);setResult("X",k,0);setResult("Y",k,0);}
	delta_avg = 1;
	delta_max = 1;
	for (n = 0;  n < 5 && delta_max > 0.01; n++) {
		delta_avg = 0;
		delta_max = 0;
		num = 0;
		for (i =  1; i <= frames-1; i++) {
			for (j = i+1; j <= minOf(i+2, frames); j++) {					
				selectImage(id);
				Stack.setFrame(i);
				run("Duplicate...", "title=A use");
				id1 = getImageID();
				selectImage(id);
				Stack.setFrame(j);
				run("Duplicate...", "title=B use");
				id2 = getImageID();
				shift = estimateShift(id1, id2);
				veloc = sqrt(shift[0]*shift[0]+shift[1]*shift[1]);
				delta_avg += veloc;
				delta_max = maxOf(veloc, delta_max);
				num++;
				//print("shift "+i+"-"+j+": "+shift[0]+", "+shift[1]+" vel:"+veloc);
				setResult("Frame", j-1, j);
				setResult("X", j-1, getResult("X", j-1)+shift[0]);
				setResult("Y", j-1, getResult("Y", j-1)+shift[1]);			
				updateResults();
				selectImage(id);
				Stack.setPosition(1, 1, j);
				run("Translate...", "x="+-shift[0]+" y="+-shift[1]+" interpolation=Bicubic slice");
				selectImage(id1);close();
				selectImage(id2);close();
				
			}					
		}
		delta_avg /= num;
		print("Residual displacement average:" + delta_avg + ", maximum:" + delta_max);		
		smoothXYResults();		
	}
	close();
}


function smoothXYResults() {
	// Smooth XY in place
	for (i = 1; i < nResults; i++) {
		setResult("X", i, (getResult("X", i-1) + getResult("X", i))/2);
		setResult("Y", i, (getResult("Y", i-1) + getResult("Y", i))/2);
		updateResults();
	}
	for (i = nResults-2; i >= 0; i--) {				
		setResult("X", i, (getResult("X", i+1) + getResult("X", i))/2);
		setResult("Y", i, (getResult("Y", i+1) + getResult("Y", i))/2);
		updateResults();
	}	
}

function estimateShift(id1,id2) {
	sx = 0;
	sy = 0;
	n = 1;
	for (j = 0; j < 2; j++) {
		scale = pow(2,j-2);
		for (iter = 0; iter < 20 && n > 0.1; iter++) {
			selectImage(id2); run("Duplicate...", " "); id3 = getImageID();
			run("Translate...", "x="+-sx+" y="+-sy+" interpolation=Bilinear");
			shift = estimateShiftCore(id1,id3,scale);
			selectImage(id3); close();
			sx -= shift[0];
			sy -= shift[1];
			n = sqrt(shift[0]*shift[0]+shift[1]*shift[1]);
		}
	}
	return newArray(sx,sy);
}

function estimateShiftCore(id1,id2,scale) {
	selectImage(id1);
	w = round(getWidth()*scale);
	h = round(getHeight()*scale);
	run("Duplicate...", " ");
	run("32-bit");
	if (scale != 1) {
		run("Size...", "width="+w+" height="+h+" depth=1 constrain average interpolation=Bilinear");	
	}
	ix = getImageID();rename("ix");
	selectImage(ix);
	run("Duplicate...", " ");
	iy = getImageID();rename("iy");

	imageCalculator("Subtract create 32-bit", id2, id1);
	it = getImageID();rename("it");
	if (scale != 1) {
		run("Size...", "width="+w+" height="+h+" depth=1 constrain average interpolation=Bilinear");	
	}

	selectImage(ix);
	run("Convolve...", "text1=[0 0 0\n-1 0 1\n0 0 0\n]");
	
	selectImage(iy);
	run("Convolve...", "text1=[0 -1 0\n0 0 0\n0 1 0\n]");
	// compute the second order images ixx iyy 
	imageCalculator("Multiply create 32-bit", ix, it); ixt = getImageID();rename("ixt");
	imageCalculator("Multiply create 32-bit", iy, it); iyt = getImageID();rename("iyt");
	imageCalculator("Multiply create 32-bit", ix, iy); ixy = getImageID();rename("ixy");
	selectImage(ix); run("Square");
	selectImage(iy); run("Square");
	// create a mask
	imageCalculator("Multiply create 32-bit", id1, id2); mask = getImageID();rename("mask");	
	getStatistics(area, mean, min, max, std, histogram);
	setThreshold(mean, max);
	run("Create Selection");
	//run("Macro...", "code=v=(v>"+mean+") slice");
	ids = newArray(ix,iy,it,ixy,ixt,iyt);
	//for (i = 0; i < ids.length; i++) {imageCalculator("Multiply", ids[i], mask);}
	
	// measure the mean intensity in each image
	getDimensions(w, h, channels, slices, frames);
	band = 2;
	selectImage(ix);  makeRectangle(band, band, w-2*band, h-2*band); run("Restore Selection"); a = getValue("Mean");
	selectImage(ixy); makeRectangle(band, band, w-2*band, h-2*band); run("Restore Selection"); b = getValue("Mean");
	selectImage(iy);  makeRectangle(band, band, w-2*band, h-2*band); run("Restore Selection"); c = getValue("Mean");
	selectImage(ixt); makeRectangle(band, band, w-2*band, h-2*band); run("Restore Selection"); d = getValue("Mean");
	selectImage(iyt); makeRectangle(band, band, w-2*band, h-2*band); run("Restore Selection"); e = getValue("Mean");
	// solve [a,b;b,c] [vx;vy] = [d;e] 
	det = a*c-b*b;
	if (abs(det) > 1e-3) {
		vx = 2*(c*d-b*e) / (a*c-b*b);
		vy = 2*(a*e-b*d) / (a*c-b*b);
	} else {
		vx = 0;
		vy = 0;
	}
	if (isNaN(vx)) {vx=0;}
	if (isNaN(vy)) {vy=0;}
	// clear the list of images
	for (i = 0; i < ids.length; i++) {selectImage(ids[i]);close();}
	//selectImage(mask);close();
	return newArray(vx/scale,vy/scale);
}


function createTestImage() {
	T = 100;
	W = 200;
	H = 200;
	d = 5;
	N = 50;	
	newImage("Test Image", "8-bit black", W, H, 1, 1, T);
	dx = 0;
	dy = 0;
	vx = 0;
	vy = 0;
	if (isOpen("ROI Manager")) { selectWindow("ROI Manager"); run("Close"); }
	setColor(255);
	bxmin = W/2+dx-d/2;
	bxmax = W/2+dx-d/2;
	bymin = H/2+dy-d/2;
	bymax = H/2+dy-d/2;
	for (t = 1; t <= T; t++) {
		setSlice(t);			
		makeOval(W/2+dx-d/2,H/2+dy-d/2, d, d);	
		roiManager("add");
		setResult("X", t-1, dx);
		setResult("Y", t-1, dy);
		setResult("Real X", t-1, dx);
		setResult("Real Y", t-1, dy);
		fill();
		bxmin = minOf(bxmin, W/2+dx-d/2);
		bxmin = minOf(bxmin, W/2+dx-d/2);
		bxmin = minOf(bxmin, W/2+dx-d/2);
		bxmin = minOf(bxmin, W/2+dx-d/2);
		for (n = 0; n < N; n++) {
			makeOval(W/2+dx-d/2+(n/5+W/4)*cos(5.5*2*PI*n/N+t),H/2+dy-d/2+(n/5+H/4)*sin(5.5*2*PI*n/N+t), d, d);		
			fill();		
		}
		dx += round(vx + 5*(random-0.5));
		dy += round(vy + 5*(random-0.5));
		vx = 0.1*random();
		vx = 0.1*random();
	}
	setSlice(1);
	makeOval(W/2-4*d,H/2-4*d, 8*d, 8*d);
}
