#@String(label="Channel names", description="list the name of the channels in the channel order, separated by comas", value="caspase,bodipy") channel_names_csv
#@String(label="ROI names", description="list the name of the ROI, separated by comas, the order will define the ROI group index", value="cell,foci") roi_names_csv
#@String(label="ROI colors", description="list the colors for the ROI, separated by comas", value="red,green") roi_colors_csv
#@String(label="Measurment", choices={"RawIntDen","Mean","Area"}, style="list") measurement

/*
 * ROI Morphing and intensity measurment
 * 
 * From a set of ROI with name defining their groups interpolate over time the ROI
 * and measure the sum of the intensity over time.
 * 
 * Jerome Boulanger 2022 for Yangci
 */

start_time = getTime();
print("\\Clear");


if (nImages==0) {
	createTest();
	testmode = true;
} else {
	testmode=false;
}

channel_names = parseCSVString(channel_names_csv);
roi_names = parseCSVString(roi_names_csv);
roi_colors = parseCSVString(roi_colors_csv);

tbl = File.getNameWithoutExtension(getTitle())+'.csv';
initTable(tbl, roi_names, channel_names);

Overlay.remove();
id = getImageID();

preprocessROIs(roi_names, roi_colors);

bloombloom(tbl, roi_names, channel_names, measurement);

addOverlays();

if (testmode) {
	createMontage();
}

end_time = getTime();
print("Elapsed time " + (end_time - start_time)/1000 + "s");

function createMontage() {
	Stack.getDimensions(width, height, channels, slices, frames);
	setColor("black");
	setFont("Arial",40);
	for (i = 1; i <= frames; i++) {		
		Overlay.drawString("t="+i, 2, 40);
		Overlay.setPosition(0, 0, i);
		Overlay.add();
	}
	Stack.setDisplayMode("composite");
	run("RGB Color", "frames keep");
	wait(100);
	id1=getImageID();
	run("Flatten", "stack");
	n = Math.ceil(sqrt(frames));
	run("Make Montage...", "columns="+n+" rows="+n+" scale=1");
	wait(100);
	id3=getImageID();
	selectImage(id1); close();	
	selectImage(id3); 
}

function bloombloom(tbl, roi_names, channel_names, measurement) {	
	// morph roi and measure over time
	for (i = 0; i < roi_names.length; i++) {
		rois0 = getROIGroup(roi_names[i]);
		rois1 = addIntermediateROIsForGroup(rois0);	
		measure(tbl, rois1, roi_names[i], channel_names,measurement);
	}
}

function addOverlays() {
	// add overlays
	n = roiManager("count");
	for (i = 0; i < n; i++) {
		roiManager("select", i);
		Roi.getPosition(channel, slice, frame);		
		Overlay.addSelection();	
		Overlay.setPosition(0, 0, frame);
		Overlay.show();
	}
}

function getROINameList() {
	// get the names of the ROI
	n = roiManager("count");
	roiManager("select", 0);
	names = newArray(Roi.getName);	
	k = 1;
	for (i = 1; i < n; i++) {
		roiManager("select", i);
		name = Roi.getName;
		for (j = 0; j < names.length; j++) {
			if (!matches(name, names[j])) {			
				names[k] = name;
				k++;
				break;
			}
		}		
	}
	return Array.trim(names, k);
}

function measure(tbl, rois, roi_name, channel_names, measurement) {
	// measure the sum of the intensity in each ROI/channel for each time point		
	for (i = 0; i < rois.length; i++) {						
		id = rois[i];
		frame = getROIFrame(id);
		for (j = 0; j < channel_names.length; j++) {			
			roiManager("select", id);	
			Stack.setChannel(j+1);	
			value = getValue(measurement);
			Table.set(roi_name+"-"+channel_names[j], frame-1, value);
		}
	}
	Table.update();
}

function initTable(tbl, roi_names, channel_names) {
	// init a table a fill it with -1;	
	Table.create(tbl);	
	n = roiManager("count");
	frames = 0;
	for (i = 0; i < n; i++) {		
		frame = getROIFrame(i);
		if (frame > frames) {
			frames = frame;
		}
	}	
	tunits = getTimeUnits();
	for (frame = 1; frame <= frames; frame++) {
		Table.set("Frame", frame-1, frame);
		Table.set(tunits[2], frame-1, (frame-1) * tunits[0]);
		for (i = 0; i < roi_names.length; i++) {
			for (j = 0; j < channel_names.length; j++) {
				Table.set(roi_names[i]+"-"+channel_names[j], frame-1, NaN);
			}
		}		
	}
}

function getTimeUnits() {
	// return an array with [ FrameInterval,Unit,"Time ["+Unit"]" ]
	dt = Stack.getFrameInterval();
	if (dt==0) {
		dt = 1;
		Time = "Frame";
	} else {
		Stack.getUnits(X, Y, Z, Time, Value);
	}
	return newArray(dt, Time, "Time ["+Time+"]");
}

function addIntermediateROIsForGroup(rois) {
	// add intermediate ROIs for a given group defined by a ordered list of 
	// roi indices.
	nrois = rois;
	for (i = 0; i < rois.length-1; i++) {
		ids = addIntermediateROIs(rois[i], rois[i+1]);				
		nrois = Array.concat(nrois, ids);		
	}
	nrois = orderROIIndicesByFrame(nrois);
	return nrois;
}

function addIntermediateROIs(id1, id2) {
	// add all missing rois between id1 and id2 based on their frame position
	t1 = getROIFrame(id1);
	name = Roi.getName;
	color = Roi.getStrokeColor;
	p1 = getROICoordinates(id1);	
	t2 = getROIFrame(id2);
	p2 = getROICoordinates(id2);
	if (t2 - t1 < 1) {return newArray();}
	D = pairwiseDistance(p1, p2);
	P = sinkhorn_uniform(D,100);	
	threshold = 0.5 * matmax(P,-1);
	print("Interpolation of roi:" +id1+"[t:"+t1+"] and roi:"+id2+" [t:"+t2+"]");
	rois = newArray(t2-t1-1);
	for (k = 0; k < rois.length; k++) {
		p = interpolateROICoordinates(P,p1,p2,t1,t2,t1+k+1,threshold);
		rois[k] = addROIfromCoordinates(p,t1+k+1,name,color);		
	}
	return rois;
}

function addROIfromCoordinates(p,t,name,color) {
	// add a oi to the roi manager from the coordinates
	n = matcol(p);
	x = newArray(n);
	y = newArray(n);
	for (i = 0; i < n; i++) {
		x[i] = matget(p,0,i);
		y[i] = matget(p,1,i);
	}	
	makeSelection(3, x, y);		
	run("Fit Spline");
	run("Interpolate", "interval=1 smooth");
	Roi.setName(name);	
	Roi.setStrokeColor(color);
	Roi.setPosition(0, 0, t);
	roiManager("add");
	id = roiManager("count")-1;
	return id;
}

function interpolateROICoordinates(P,p0,p1,t0,t1,t,threshold) {
	// interpolate ROI coordinates between two time points t0 and t1 at time t
	n = matcol(p0);
	m = matcol(p1);
	p = matzeros(2,n*m);
	dt = (t-t0) / (t1-t0-1);
	// this make the deformation more progressive at the begining and the end
	s = 3 * dt * dt - 2 * dt * dt * dt;		
	k = 0;
	for (i = 0; i < n; i++) {
		for (j = 0; j < m; j++) {			
			if (matget(P,i,j) > threshold) {
				matset(p,0,k,(1-s) * matget(p0,0,i) + s * matget(p1,0,j));
				matset(p,1,k,(1-s) * matget(p0,1,i) + s * matget(p1,1,j));			
				k = k + 1;
			}
		}
	}
	// crops the coordinate columns to [0-k]
	pfinal = matzeros(2,k-1);
	for (i = 0; i < k-1; i++) {
		matset(pfinal,0,i,matget(p,0,i));
		matset(pfinal,1,i,matget(p,1,i));
	}
	return pfinal;
}


function orderROIIndicesByFrame(rois) {
	// order the ROIs indices in rois by their frame number
	time = newArray(rois.length);
	for (i = 0; i < rois.length; i++) {
		roiManager("select", rois[i]);
		wait(10);
		Roi.getPosition(channel, slice, frame);
		time[i] = frame;		
	}		
	ranks = Array.rankPositions(time);
	ordered_rois = newArray(rois.length);
	for (i = 0; i < rois.length; i++) {
		ordered_rois[i] = rois[ranks[i]];
	}	
	return ordered_rois;
}

function getROIFrame(i) {
	// Frame position of the ROI at index id
	roiManager("select", i);
	Roi.getPosition(channel, slice, frame);
	return frame;
}

function getROIGroup(name) {
	// list the roi from one group using their name and sort them by frame index
	n = roiManager("count");
	rois = newArray(n);
	time = newArray(n);
	k = 0;	
	for (i = 0; i < n; i++) {
		roiManager("select", i);		
		if (matches(Roi.getName(), name)) {
			Roi.getPosition(channel, slice, frame);
			rois[k] = i;
			time[k] = frame;
			k++;
		}
	}
	rois = Array.trim(rois, k);
	time = Array.trim(time, k);
	ranks = Array.rankPositions(time);
	ordered_rois = newArray(rois.length);
	for (i = 0; i < rois.length; i++) {
		ordered_rois[i] = rois[ranks[i]];
	}
	return ordered_rois;
}

function preprocessROIs(names,colors) {
	// assign groups to the ROI and interpolate the individual contours
	for (i = 0; i < names.length; i++) {
		names[i] = toLowerCase(names[i]);
	}	
	// lowercase name, assign groups, colors and interpolate
	n = roiManager("count");
	for (i = 0; i < n; i++) {
		roiManager("select", i);
		roiManager("rename", String.trim(toLowerCase(Roi.getName())));
		for (k = 0; k < names.length; k++) {
			if (matches(Roi.getName(), names[k])) {
				RoiManager.setGroup(k+1);
				Roi.setStrokeColor(colors[k]);
				roiManager("update");
				break;
			}
		}
		getPixelSize(unit, px, px);
		ds = maxOf(2, getValue("Perim.") / px / 50);		
		print(ds);
		run("Interpolate", "interval="+ds+" smooth");
		run("Fit Spline");
		run("Interpolate", "interval="+ds+" smooth");
		roiManager("update");
	}
}

function createTest() {
	print("Generating a test image");
	run("Roi Defaults...", "color=green stroke=1 group=1 name=cell");
	run("Roi Defaults...", "color=red stroke=1 group=2 name=foci");
	if (isOpen("ROI Manager")) {selectWindow("ROI Manager");run("Close");}
	width = 400;
	height = 400;
	channels = 3;
	frames = 16;
	newImage("Test Image", "8-bit white", width, height, channels, 1, 2);
	Stack.setDisplayMode("composite");
	Stack.setPosition(1, 1, 1);
	makeRectangle(width/4, height/4, width/2, height/2);
	setColor(255);
	run("Fill", "slice");
	Stack.setPosition(1, 1, 2);
	makeOval(width/4, height/4, width/2, height/2);
	setColor(128);
	run("Fill", "slice");
	Stack.setPosition(2, 1, 1);
	makeOval(0.6*width, 0.6*height, 0.1*width, 0.1*width);
	setColor(255);
	run("Fill", "slice");
	Stack.setPosition(2, 1, 2);
	makeOval(0.8*getWidth, 0.8*getHeight, 0.12*getWidth, 0.12*getWidth);	
	setColor(255);
	
	run("Fill", "slice");
	run("Select None");	
	run("Size...", "width="+width+" height="+height+" time="+frames+" average interpolation=Bilinear");
	
	// define cells at 3 time points
	Stack.setFrame(1);
	makeRectangle(width/4, height/4, width/2, height/2);	
	run("Interpolate", "interval=2 smooth");
	Roi.setGroup(1);
	Roi.setName("cell");
	Roi.setPosition(1, 1, 1);
	roiManager("add");
		
	Stack.setFrame(frames);
	makeOval(width/4, height/4, width/2, height/2);
	run("Interpolate", "interval=2 smooth");	
	Roi.setGroup(1);
	Roi.setName("cell");
	Roi.setPosition(1, 1, frames);
	roiManager("add");	
	
	Stack.setFrame(round(frames/2));	
	// draw a star
	xpts = newArray(10);
	ypts = newArray(10);
	for (k = 0; k < 10; k++) {
		xpts[k] = getWidth*(0.5+(0.4*(k%2==0)+0.2*((k+1)%2==0))*cos(2*PI*k/10));
		ypts[k] = getHeight*(0.5+(0.4*(k%2==0)+0.2*((k+1)%2==0))*sin(2*PI*k/10));
	}
	makeSelection("polygon",xpts,ypts);
	run("Interpolate", "interval=2 smooth");
	Roi.setGroup(1);
	Roi.setName("cell");
	Roi.setPosition(1, 1, round(frames/2));
	roiManager("add");
		
	// define foci at 2 time points
	Stack.setFrame(1);
	makeOval(0.4*width, 0.4*height, 0.1*width, 0.1*width);	
	Roi.setGroup(2);
	Roi.setName("foci");
	Roi.setPosition(2, 1, 1);
	roiManager("add");
	
	Stack.setFrame(frames);		
	makeOval(0.6*getWidth, 0.6*getHeight, 0.12*getWidth, 0.12*getWidth);	
	Roi.setGroup(2);
	Roi.setName("foci");
	Roi.setPosition(2, 1, frames);
	roiManager("add");	
	channel_names_csv = "caspase,bodipy";
	roi_names_csv = "cell,foci";
}

function parseCSVString(csv) {
	// split a csv string and trim each element
	str = split(csv, ",");
	for (i = 0; i < str.length; i++) { 
		str[i] =  String.trim(str[i]);
	}
	return str;
}


function diff(u) {
	// centered finite difference with symmetric bc
	n = u.length;
	d = newArray(n);	
	d[0] = u[1] - u[0];	
	for (i = 1; i < n - 1; i++) {
		d[i] = (u[i+1] - u[i-1]) / 2.0;
	}
	d[n-1] = u[n-1] - u[n-2];
	return d; 
}


function sinkhorn_uniform(C,gamma) {
	// Optimal transport using uniform weights
	n = matrow(C);
	m = matcol(C);
	p = matones(n,1);
	q = matones(m,1);
	P = sinkhorn(C,p,q,gamma);
	return P;
}

function sinkhorn(C,p,q,gamma) {
	/* Optimal transport 
	 *  Parameters
	 * X: cost [n,m] matrix (distance)
	 * p : mass [m,1] vector (probability)
	 * q : mass [m,1] vector (probability)
	 * gamma : entropic regularization parameter
	 * Note
	 * - p and q will be normalized by their sum and reshaped as column vector
	 * - if gamma < 0, gamma is set to 10/min(C) 
	 * 
	 */
	n = matrow(C);
	m = matcol(C);
	assert((matsize(p) == n), "Size of C ["+n+"x"+m+"] and p ["+matsize(p)+"] do not match");
	assert((matsize(q) == m), "Size of C ["+n+"x"+m+"] and p ["+matsize(q)+"] do not match");
	p = matreshape(p,n,1);
	q = matreshape(q,m,1);
	// normalize p and q
	p = mat_normalize_sum(p);
	q = mat_normalize_sum(q);
	
	// compute a default gamma value from C
	if (gamma < 0) {
		minc = C[2];
		for (i = 3; i < C.length; i++) {minc = minOf(minc, C[i]);}		
		gamma = 10*minc;
		print("min c: "+minc+" gamma = "+gamma);
	}
	// compute xi = exp(-C/gamma)
	xi = matzeros(n,m);	
	for (i = 0; i < n; i++) {
		for (j = 0; j < m; j++) {
			cij = matget(C,i,j);
			val = exp(-cij/gamma);
			matset(xi,i,j,val);
		}
	}
	// iterativerly evaluate
	// a = p ./ (xi*b)
	// b = q ./ (x'*a)
	xit = mattranspose(xi);
	b = matones(m,1);
	for (iter = 0; iter < 30; iter++) {
		xib = matmul(xi,b);
		a = mat_div(p,xib);
		xita = matmul(xit,a);
		b = mat_div(q,xita);
	}
	// P = diag(a)*xi*diag(b)
	a = matdiag(a);
	b = matdiag(b);
	P = matmul(xi, b);
	P = matmul(a, P);
	
	return P;
}

function mat_normalize_sum(A) {
	// normalize A by the sum over all indices
	B = A;
	s = 0;
	for (i = 2; i < A.length; i++) {
		 s += A[i];
	}
	for (i = 2; i < A.length; i++) {
		 B[i] = A[i] / s;
	}
	return B;
}

function assert(cond, msg) {
	// check that condition is true
	if (!(cond)) {
		exit(msg);
	}
}

function pairwiseDistance(p,q) {
	// Compute the pairwie distance between p and q
	// P : 2xn coordinate matrix
	// Q : 2xm coordinate matrix
	// return distance [n x m] matrix
	n = matcol(p);
	m = matcol(q);
	D = matrow(p);
	a = matzeros(n,m);
	for (i = 0; i < n; i++) {
		for (j = 0; j < m; j++) {
			s = 0;
			for (k = 0; k < D; k++) {
				d = matget(p,k,i) - matget(q,k,j);
				s += d*d;
			}
			matset(a,i,j,s);
		}
	}
	return a;
}


function getROICoordinates(roi_id) {
	// get the spline anchor of the ROI
	// return a 2xn matrixs
	roiManager("select",roi_id);	
	Roi.getSplineAnchors(x, y);
	n = x.length;
	vals = Array.concat(x,y);
	m = matnew(2,n,vals);
	return m;
}


/////////////////////////////////////////////////////////
// Basic Linear Algebra
////////////////////////////////////////////////////////
function matnew(n,m,vals) {
	// Initialize a matrix with values
	// Column major, N rows and M columns
	// Matrices are array with 2 first elements being the dimensions
	A = newArray(2+n*m);
	A[0] = n;
	A[1] = m;
	for (i = 0; i < vals.length; i++) {
		A[2+i] = vals[i];
	}
	return A;
}

function matrow(A) {
	// return the number of rows
	return A[0];
}

function matcol(A) {
	// return the number of columns
	return A[1];
}

function matsize(A) {
	// return the number of elements of A
	return A.length - 2;
}

function matzeros(n,m) {
	// Return an NxM matrix initialized with zeros
	A = newArray(2+n*m);
	A[0] = n;
	A[1] = m;
	return A;
}

function matones(n,m) {
	// Return an NxM matrix initialized with zeros
	A = matzeros(n,m);
	for (i = 2 ; i < A.length; i++) {
		A[i] = 1;
	}
	return A;
}

function matreshape(A,n,m) {
	// Reshape the matrix
	B = A;
	B[0] = n;
	B[1] = m;
	return B;
}

function matdiag(v) {
	// create a diagonal matrix with the values of the vector v
	n = v.length - 2;
	A  = matzeros(n,n);
	for (i = 0 ; i < n; i++) {
		matset(A,i,i,matget(v,i,0));
	}
	return A;
}

function matset(A,i,j,val) {
	// set value of A[i,j] in place
	if (i >= 0 && i < A[0] && j >=0 && j < A[1]) {
		A[2+j+A[1]*i] = val;
	} else {
		print("Error in matset: out of bound [" + i + ","+j+"] in A["+A[0]+","+A[1]+"]");
		please_check_the_message_in_the_log();
	}
}

function matget(A,i,j) {
	// get the value A[i,j]
	if (i >= 0 && i < A[0] && j >=0 && j < A[1]) {
		return A[2+j+A[1]*i];
	} else {
		print("Error in matset: out of bound [" + i + ","+j+"] in A["+A[0]+","+A[1]+"]");
		please_check_the_message_in_the_log();
	}
}

function mattranspose(A) {
	// transpose matrix A
	M = matcol(A);
	N = matrow(A);
	B = matzeros(M,N);
	for (n=0;n<N;n++) {
		for (m=0;m<M;m++) {
			B[2+n+m*N] = A[2+m+n*M];
		}
	}
	return B;
}

function mat2str(A,name) {
	// Convert the matrix to a string
	if (A.length == 0){
		str = name + "is an empty matrix";
	} else {
		N = A[0];
		M = A[1];
		str = name + " ["+N+"x"+M+"] :\n";
		for (n=0;n<N;n++) {
			for (m=0;m<M;m++) {
				str = str + A[2+m+n*M];
				if ( m == M-1 ){ str = str + "\n";}
				else  {str = str + ", ";}
			}
		}
	}
	return str;
}

function mat2img(A,name) {
	// create an image and fill it with the matrix values
	width = matcol(A);
	height = matrow(A);
	newImage(name, "32-bit", width, height, 1);
	id = getImageID();
	for (y = 0; y < height; y++) {
		for (x = 0; x < width; x++) {
			setPixel(x, y, matget(A,y,x));
		}
	}
	run("Enhance Contrast", "saturated=0.0");
	run("In [+]");
	run("In [+]");
	return id;
}

function matadd(A,B,alpha) {
	// Add 2 matrices A and B (A + alpha B)
	C = matzeros(A[0],A[1]);
	for (i = 2; i <A.length; i++){C[i] = A[i]+alpha*B[i];}
	return C;
}

function matrep(A,K,L) {
	// repeat matrix A
	N = A[0];
	M = A[1];
	B = matzeros(N*K,L*M);
	for (k = 0; k < K; k++) {
		for (l = 0; l < L; l++) {	
			for (n = 0; n < N; n++) {
				for (m = 0; m < M; m++) {
					B[2+m+l*M+(M*L)*(n+k*N)] = A[2+m+n*M];
				}
			}
		}
	}
	return B;
}

function matmul(A,B) {
	// matrix multiplication C[nxk] = A[nxm] x B[mxk] 
	// C[n,k] = sum A[n,m] B[m,k]
	// (column major matrices)
	N = A[0];
	M = A[1];
	K = B[1];
	if (B[0] != A[1]) {
		print("matmul(A,B): incompatible dimensions A ["+ A[0]+"x"+A[1]+"], B ["+B[0]+"x"+B[1]+"].");
		please_check_the_message_in_the_log();
		
	}
	C = matzeros(N,K);
	for (n = 0; n < N; n++) {
		for (k = 0; k < K; k++) {
			for (m = 0; m < M; m++) {
				C[2+k+n*K] = C[2+k+n*K] + A[2+m+n*M] * B[2+k+m*K];
			}
		}
	}
	return C;
}

function mat_mul(A,B) {
	// multiply entry wise the two matrices
	// C[i,j] = A[i,j] * B[i,j]
	assert((A.length == B.length), "A and B have incompatible size");
	N = A[0];
	M = A[1];
	C = matzeros(N,M)
	for (i = 2; i < C.length; i++) {
		C[i] = A[i] * B[i];
	}
	return C;
}

function mat_div(A,B) {
	// divide entry wise A by B
	// C[i,j] = A[i,j] / B[i,j]
	assert((A.length == B.length), "A and B have incompatible size");
	N = A[0];
	M = A[1];
	C = matzeros(N,M);
	for (i = 2; i < C.length; i++) {
		C[i] = A[i] / B[i];
	}
	return C;
}


function matsum(A,dim) {
	// sum values of A algon 1 dimension (0 or 1), if dim==-1, sum over all values
	N = A[0];
	M = A[1];
	if(dim==0) {
		S = matzeros(1,M);
		for (n = 0; n < N; n++) {
			for (m = 0; m < M; m++) {
				S[2+m] = S[2+m] + A[2+m+n*M];
			}
		}
	} else if (dim==1){
		S = matzeros(N,1);
		for (n = 0; n < N; n++) {
			for (m = 0; m < M; m++) {
				S[2+n] = S[2+n] + A[2+m+n*M];
			}
		}
	} else {
		S = 0;
		for (n = 0; n < N; n++) {
			for (m = 0; m < M; m++) {
				S = S + A[2+m+n*M];	
			}
		}
	}
	return S;
}

function matmean(A,dim) {
	// average along dimension if dim==0 or 1
	B = matsum(A,dim);
	if (dim == 0 || dim==1) {
		for (i = 2; i < B.length; i++) {
			B[i] = B[i] / A[dim];
		}
	} else {
		B = B / (A[0]*A[1]);
	}
	return B;
}

function matmax(A,dim) {
	// average along one dimension
	N = A[0];
	M = A[1];
	if(dim==0) {
		S = matzeros(1,M);
		for (n = 0; n < N; n++) {
			for (m = 0; m < M; m++) {
				S[2+m] = maxOf(S[2+m],A[2+m+n*M]);
			}
		}
	} else if (dim==1){
		S = matzeros(N,1);
		for (n = 0; n < N; n++) {
			for (m = 0; m < M; m++) {
				S[2+n] = maxOf(S[2+n],A[2+m+n*M]);
			}
		}
	} else {
		S = A[2];
		for (n = 0; n < N; n++) {
			for (m = 0; m < M; m++) {
				S = maxOf(S,A[2+m+n*M]);	
			}
		}
	}
	return S;
}

function solve(A,B) {
	// Solve A*X=B by iterating X=X-At(AX - B)
	// A[NxM] x X[MxK] = B [NxK]
	print(mat2str(A,"A"));
	print(mat2str(B,"B"));
	N = A[0];
	M = A[1];
	K = B[1];
	At = mattranspose(A);
	AtA = matmul(At,A);
	AtB = matmul(At,B);
	X = matzeros(M,K);
	a = matnorm(AtA);
	for (iter = 0; iter < 1000; iter++) {
		AtAX = matmul(AtA,X);
		dX = matadd(AtAX,AtB,-1);
		X = matadd(X, dX, -1/(2*a));
	}
	return X;
}

function matnorm(A) {
	// L2 norm of the matrix (sum a_ij^2)
	n = 0;
	for (i = 2; i < A.length; i++) {
		n = n + A[i]*A[i];
	}
	return sqrt(n);
}

