/*
 * FRAP intensity measure with manual motion correction
 * 
 * Add several regions of interest to the ROI manager:
 *  - create a roi, optimize its position in x,y and pressing the keys g and z
 *  - add the ROI to the roi maanger (pressing u)
 *  
 *  Run the macro (press q)
 *  - set the radius for grouping ROIs
 *  - set the diameter of the measurement region
 *  
 *  Load the image with bioformat to get the correct frame interval froma LSM file.
 *  
 * Jerome Boulanger 2020 for Yangci
 */


macro "Measure Intensities [q]" {	
	Dialog.create("FRAP Measure intensities");
	Dialog.addNumber("Grouping radius [px]", 50);
	Dialog.addNumber("Measurement diameter [px]", 10);
	Dialog.show();
	R = Dialog.getNumber();
	region_size = Dialog.getNumber();	
	cmap = newArray("#aa0000","#00aaaa","#00aa00","#aa00aa","#0000aa","#aabb55");	
	id = getImageID;
	Overlay.remove;
	setBatchMode(true);
	data = trackROI(R);		
	analyseFrapRecovery(data);
	setBatchMode(false);
	selectImage(id);
}


macro "Optimize XY position [g]" {
	print("\\Clear");		
	getVoxelSize(dx, dy, dz, unit);	
	for (i = 0 ; i <10; i++) {
		Roi.getBounds(x, y, width, height);		
		xc = getValue("X");
		yc = getValue("Y");
		xm = getValue("XM");
		ym = getValue("YM");
		vx = (xm-xc)/dx;
		vy = (ym-yc)/dy;
		Roi.move(x+vx,y+vy);
	}	
	//roiManager("update");
}


macro "Find Brightest Z plane [z]" {
	Stack.getDimensions(width, height, channels, slices, frames);
	imax = 0;
	zbest = 1;
	for (n = 1; n < slices; n++) {
		Stack.setSlice(n);
		v = getValue("Mean");
		print("z="+n+" v="+v);
		if (v > imax) {
			zbest = n;
			imax = v;
		}
	}
	Stack.setSlice(zbest);
	//roiManager("update");
}


function analyseFrapRecovery(data) {
	Overlay.remove;
	run("Select None");
	n = getNumEntry(data);
	L = getNumLabels(data);
	Stack.getUnits(X, Y, Z, Time, Value);
	dt = Stack.getFrameInterval();
	if (dt == 0) {
		print("Uncalibrated time interval");
		dt = 1;
	}
	// Display the ROI labels on the stack
	labels_str = newArray(L);	
	for ( i = 0; i < L; i++) {
		labels_str[i] = "ROI "  + (i + 1);
		for (j = 0; j < n; j++) {
			if (getVal(data,j,"t")==1 && getVal(data,j,"label")==i) {				
				Overlay.drawString(labels_str[i], getVal(data,j,"x"), getVal(data,j,"y") );	
				Overlay.setPosition(1,0,0);								
			}
		}		
	}
	Overlay.show();
	
	//  ask for identification of the ROI
	Stack.setPosition(1,0,0);
	txt = newArray(L);
	perm = newArray(L);
	for (k = 0; k < L; k++) {		
		if (k == 0) { 
			txt[k] = "FRAP";
		} else if (k == 1) {
			txt[k] = "Background";
		} else if (k == 2) {
			txt[k] = "Control";
		} else {
			txt[k] == "Other";
		}
	}
	Dialog.create("Frap recovery:ROI selection");
	Dialog.addNumber("FRAP frame", 3);	
	for (k = 0; k < L; k++) {	
		Dialog.addChoice(txt[k], labels_str, labels_str[k]);
	}	
	Dialog.show();
	frap_frame = Dialog.getNumber();
	Overlay.remove();
	for (i = 0; i < L; i++) {	
		str = Dialog.getChoice();	
		for (j = 0; j < L; j++) {
			if (matches(labels_str[j], str)) {				
				perm[i] = j;
				for (k = 0; k < n; k++) {				
					if (getVal(data,k,"t")==1 && getVal(data,k,"label")==j) {						
						setColor(cmap[i]);
						Overlay.drawString(txt[i], getVal(data,k,"x"), getVal(data,k,"y") );	
						Overlay.setPosition(1,0,0);	
					}
				}
			}
		}		
	}
	Overlay.show();	
	
	// export data to a table
	createFrapResultTable(data,perm,txt);
	
	// create a graph
	Plot.create("Intensity measure", "Time ["+Time+"]", "Intensity");
	Plot.setLineWidth(5);
	for (l = 0; l < perm.length; l++) {		
		Plot.setColor(cmap[perm[l]]);
		X = getCurve(data, perm[l], "t");
		for (i = 0 ; i < X.length; i++) {
			X[i] = X[i] * dt;
		}
		Y = getCurve(data, perm[l], "i");		
		Plot.add("line", X, Y);		
	}
	Plot.setLegend("Frap\nBackground\nControl");
	Plot.setLimitsToFit();
	Plot.update();

	a = newArray(1,2,3,4,5,6,7,8,9,10);
	a = Array.slice(a,4,10);
	Array.print(a);

    // Final graph with curve fitting
	if (dt == 0) {
		print("Uncalibrated time interval");
		dt = 1;
	}
    Tm = getCurve(data, perm[0], "t");
    Im = getCurve(data, perm[0], "i");
    Bm = getCurve(data, perm[1], "i");
    Cm = getCurve(data, perm[2], "i");    
    I0 = Im[frap_frame-1];
    I1 = Im[frap_frame];
    print("Intenity before bleaching :" + I0);
    print("Intenity after bleaching :" +  Im[frap_frame]);
    T = Array.slice(Tm, frap_frame, Tm.length-1);
    I = Array.slice(Im, frap_frame, Im.length-1);
    B = Array.slice(Bm, frap_frame, Bm.length-1);
    C = Array.slice(Cm, frap_frame, Cm.length-1);
    for (i = 0; i < T.length; i++) {
    	T[i] = (T[i] - frap_frame - 1) * dt;
    	C[i] = C[i] - B[i];    	 	
    }
    Fit.doFit("Exponential",T,C);
    bleaching_amp = Fit.p(0);
	bleaching_rate = Fit.p(1);
	print("Acquisition bleaching amplitude " + bleaching_amp);
	print("Acquisition bleaching rate " + bleaching_rate);
	for (i = 0; i < T.length; i++) {
		I[i] = (I[i] - I1) / (I0 - I1) * exp(-i*dt*bleaching_rate);	
	}
	Fit.doFit("Exponential Recovery (no offset)",T,I);
	rec_amp = Fit.p(0);
	rec_rate = Fit.p(1);
	print("Recovery rate: " + rec_rate + "/"+Time);
	print("Recovery amplitude: " + rec_amp*100 + "%");

	M = newArray(T.length);
	for (i = 0; i < M.length; i++) {
		M[i] = Fit.f(T[i]);
	}
	Plot.create("FRAP recovery ","Time ["+Time+"]", "Normalized Intensity");
	Plot.setLineWidth(5);
	Plot.setColor(cmap[0]);	
	Plot.add("cricle", T, I);
	Plot.setColor(cmap[1]);	
	Plot.add("line", T, M);
	Plot.setLimitsToFit();
	Plot.setLegend("Data\nModel");
	Plot.update();

	Array.show(T,I);
}


// group ROI, interpolate the position in time and measure the intensity
// along the position over time
function trackROI(R) {
	id  = getImageID();
	print("\\Clear");	
	print("ROI grouping");
	x = roi2feat();	
	L = groupFeat(x,R);
	print(" -number of groups : " + L);
	colorROIbyLabel(x);	
	print("Interpolate and Measure");
	run("Select None");
	y = interpolate(x,L);
	return y;	
}

function createFrapResultTable(src,perm,txt) {
	if (isOpen("Results")) {
		selectWindow("Results");
		run("Close");
	}
	n = getNumEntry(src);
	dt = Stack.getFrameInterval();	
	if (dt == 0) {
		print("Uncalibrated time interval");
		dt = 1;
	}
	Stack.getUnits(X, Y, Z, Time, Value);
	t = getCurve(src,0,"t");
	for (i = 0; i < n; i++) {
		t = getVal(src,i,"t") - 1;
		I = getVal(src,i,"i");
		l = getVal(src,i,"label");		
		// find the column corresponding to this label
		k = -1;
		for (j = 0; j < perm.length; j++)	 {
			if (perm[j] == l) {
				k = j;
			}
		}
		if (k>=0) {
			setResult("T["+Time+"]", t, t*dt);
			setResult(txt[k], t, I);			
		}
	}	
}
/*
function createFrapResultTable(src) {
	if (isOpen("Results")) {
		selectWindow("Results");
		run("Close");
	}
	n = getNumEntry(src);
	dt = Stack.getFrameInterval();	
	if (dt == 0) {
		print("Uncalibrated time interval");
		dt = 1;
	}
	Stack.getUnits(X, Y, Z, Time, Value);
	t = getCurve(src,0,"t");
	for (i = 0; i < n; i++) {
		t = getVal(src,i,"t") - 1;
		I = getVal(src,i,"i");
		l = getVal(src,i,"label");		
		setResult("T["+Time+"]", t, t*dt);		
		setResult("I"+(l+1), t, I);
	}	
}
*/
// get a values over time for a specific label and field
function getCurve(src,label,field) {		
	n = getNumEntry(src);
	frames = 0;
	for (i = 0; i < n; i++) {
		frames = maxOf(frames, getVal(src,i,"t"));
	}
	dst = newArray(frames);
	ts = newArray(frames);	
	j = 0;
	for (i = 0; i < n; i++) {
		li = getVal(src,i,"label");
		//print(li);
		if (li == label) {
			//print("ok " + li +"=="+label+" "+j+"/"+frames);
			if (j < frames) {
				dst[j] = getVal(src,i,field);
				ts[j]  = getVal(src,i,"t");				
				j = j + 1;				
			}
		}
	}
	// order the array by time points	
	r = Array.rankPositions(ts);
	dsto = newArray(frames);
	for (i = 0; i < frames; i++) {
		dsto[i] = dst[r[i]];
	}
	return dsto;
}


// interpolate data in x with L labels over regular time points
// and measure intensities at the interpolated positions
function interpolate(x,L) {	
	Stack.getDimensions(width, height, channels, slices, frames);
	y = newArray(L*frames*6);
	n = getNumEntry(x);
	
	for (l = 0; l < L; l++) { // for each label
		// determine first and last frame of the ROIs
		tstart = 0;
		istart = 0; // first element of x from label l
		tend = 0;
		iend = 0; // last element of x from label l
		for (i = 0; i < n; i++) {
			li = getVal(x,i,"label");			
			if (li == L) {
				t = getVal(x,i,"t");
				if (t <= tstart) {
					tstart = t;
					istart = i;
				}
				if (t >= tend) {
					tend = t;
					iend = i;
				}
			}
		}		
		for (frame = 1; frame <= frames; frame++) { // for each frame
			// find the control points
			dt0 = -frames;
			i0 = -1;			
			dt1 = frames;
			i1 = -1;			
			for (i = 0; i < n; i++) {	
				li = getVal(x,i,"label");			
				if (li == l) {
					dt = getVal(x,i,"t") - frame;					
					if (dt <= 0 && dt > dt0) {						
						dt0 = dt;
						i0 = i;
					}
					if  (dt >= 0 && dt < dt1) {
						dt1 = dt;
						i1 = i;
					}
				}
			}					
			// interpolate
			if (i0 != -1 && i1 != -1) {
				//print("interpolating for frame " + frame);			
				xj = interp1val(getVal(x,i0,"t"), getVal(x,i0,"x"), getVal(x,i1,"t"), getVal(x,i1,"x"), frame);
				yj = interp1val(getVal(x,i0,"t"), getVal(x,i0,"y"), getVal(x,i1,"t"), getVal(x,i1,"y"), frame);
				zj = interp1val(getVal(x,i0,"t"), getVal(x,i0,"z"), getVal(x,i1,"t"), getVal(x,i1,"z"), frame);
			} else {
				if (i0 == -1) { // before the first
					xj = getVal(x,istart,"x");
					yj = getVal(x,istart,"y");
					zj = getVal(x,istart,"z");					
				}
				if (i1 == -1) { // after the last
					xj = getVal(x,iend,"x");
					yj = getVal(x,iend,"y");
					zj = getVal(x,iend,"z");
				}
			}
			j = (frame-1) + l * frames;
			setVal(y,j,"x", xj);
			setVal(y,j,"y", yj);
			setVal(y,j,"z", zj);
			setVal(y,j,"t", frame);			
			run("Select None");
			Stack.setPosition(1, round(zj), frame);
			makeOval(xj-region_size/2,yj-region_size/2,region_size,region_size);	
			setVal(y,j,"i", getValue("Mean"));
			setVal(y,j,"label", l);			
		}
	}
	return y;
}

function interp1val(x1,y1,x2,y2,xq) {
	val = 0;
	if (abs(x1-x2) < 0.5) {
		val = y1;
	} else {
		val = y1 + (xq-x1) * (y2-y1) / (x2-x1);
	}
	return val;
}

function colorROIbyLabel(x) {	
	n = getNumEntry(x);	
	for (i = 0; i < n; i++)  {
		roiManager("select",i);		
		L = getVal(x,i,"label");
		if (L >= 0) {
			Roi.setStrokeColor(cmap[L]);		
		}
	}
}

function groupFeat(x,R) {	
	print("Grouping radius is " + R);
	n = getNumEntry(x);
	L = 0;
	D = computeDistanceMatrix(x);		
	setVal(x,0,"label",0);
	for (i = 0; i < n; i++)  {
		li = getVal(x,i,"label");		
		dmin = 1e6;
		jmin = -1; 			
		for (j = 0; j < n; j++) if (i!=j) {		
			lj = getVal(x,j,"label");
			if (D[i+n*j] < dmin && lj != -1) {
				dmin = D[i+n*j];
				jmin = j;
			}
		}
		//print(i + "-" + jmin + ", d =" + dmin);				
		if (dmin < R && jmin >= 0) {			
			lj = getVal(x,jmin,"label");
			setVal(x,i,"label",lj);			
		} else {
			setVal(x,i,"label",L);
			L = L + 1;
		}					
	}
	return L;
}

function computeDistanceMatrix(x) {
	n = getNumEntry(x);
	D = newArray(n*n);
	for (i = 0; i < n; i++)  {
		xi = getVal(x,i,"x");
		yi = getVal(x,i,"y");
		zi = getVal(x,i,"z");
		ti = getVal(x,i,"t");
		for (j = 0; j < n; j++) {
			dx = getVal(x,j,"x") - xi;							
			dy = getVal(x,j,"y") - yi;				
			dz = getVal(x,j,"z") - zi;			
			D[i+n*j] = sqrt(dx*dx+dy*dy+dz*dz);
		}
	}
	return D;
}

function printFeat(X){
	n = getNumEntry(X);
	for (i = 0; i < n; i++)  {
		setResult("x", i, getVal(X,i,"x"));
		setResult("y", i, getVal(X,i,"y"));
		setResult("z", i, getVal(X,i,"z"));
		setResult("t", i, getVal(X,i,"t"));
		setResult("i", i, getVal(X,i,"i"));
		setResult("label", i, getVal(X,i,"label"));			
	}
	updateResults();
}
// return a 2D ARRAY x,y,z,t,i,label for each roi in pixels
function roi2feat() {	
	n = roiManager("count");
	N = 6*n;
	x = newArray(N);// x,y,z,t,i,label nx6		
	getVoxelSize(dx, dy, dz, unit)
	for (i = 0; i < n; i++) {		
		roiManager("select", i);
		Stack.getPosition(channel, slice, frame);		
		List.setMeasurements;
		xi = List.getValue("X");
		yi = List.getValue("Y");
		mi = List.getValue("Mean");		
		setVal(x,i,"x",xi / dx);
		setVal(x,i,"y",yi / dy);
		setVal(x,i,"z",slice);		
		setVal(x,i,"i",mi);
		setVal(x,i,"t",frame);	
		setVal(x,i,"label",-1);	
	}
	//printFeat(x);
	return x;
}

function setVal(x,i,field,value) {
	idx = getIndex(x,i,field);	
	x[idx] = value;
	return x;
}

function getVal(x,i,field) {
	idx = getIndex(x,i,field);	
	value = x[idx];	
	return parseFloat(value);
}

function getNumEntry(x) {
	return x.length/6;
}

function getNumLabels(x) {
	lmax = 0;
	n = getNumEntry(x);
	for (i = 0; i < n; i++) {
		l = getVal(x,i,"label");
		if (l>L) {
			lmax = l;
		}
	}
	return lmax+1;
}

function getIndex(x,i,field) {
	a = newArray("x","y","z","t","i","label");
	j = 0;
	for (k = 0; k < a.length; k++) {
		if (matches(a[k], field)) {
			j = k;	
		}
	}
	idx = j+6*i;
	if (idx >= x.length || idx < 0) {print("Out of bound ("+i+","+field+")");}
	return idx;
}