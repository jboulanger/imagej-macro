#@String(choices={"Estimate Beads","Estimate ROI","Apply","Show beads"}) mode
#@String(choices={"Affine","Quadratic"}) model
#@String(label="table name",value="Models.csv") tblname

/*
 * Estimate and apply channel alignement from a bead calibration image
 * 
 * Channels are aligned to the 1st channel
 * 
 * Estimate transform:
 *  open an image with beads and save the Models.csv file
 *  
 * Apply transform:
 *   parameters are loaded from a Models.csv table
 *   
 *   
 *  Toubleshouting
 *  
 *  	1. The name of the model table has to be 'Models.csv'. Rename if needed
 *  	2. The macro will work on small images (tested on 1024x1024) and is relatively slow.
 *  	3. Some error may occur, try closing the results and models.csv tables.
 *   
 * Jerome Boulanger 2021-2022 for Sina, Anna & Ross
 */

if (nImages==0) {
	path = File.openDialog("Open an image");
	open(path);
}

if (matches(mode, "Estimate Beads")) {	
	print("Estimate transfrom from beads");
	estimateTransform(model,tblname);
} else if (matches(mode, "Estimate ROI")) {
	print("Estimate transfrom from ROIs");
	estimateTransformROI(model,tblname);
} else if (matches(mode,"Show beads")) {
	showBeads();
} else {	
	print("Apply transform");
	applyTransform(tblname);
}
print("Done");

function showBeads(){
	// load coordinates from the result table and add them to the ROI manager
	n = nResults;
	for (i = 0; i < n; i++) {
		x = getResult("X", i);
		y = getResult("Y", i);
		c = getResult("Channel", i);		
		makePoint(x, y);
		Roi.setPosition(c, 1, 1);
		roiManager("add");
	}
	roiManager("show all without labels");
}

function applyTransform(tblname) {	
	Stack.getDimensions(width, height, channels, slices, frames);
	for (channel = 2; channel <= channels; channel++) {
		cname = "M1"+channel;
		print("Loading parameters from column " + cname);
		if (isOpen(tblname)) {
			selectWindow(tblname);
		} else {
			exit("Could not find table " + tblname);			
		}
		M = Table.getColumn(cname);
		for (frame = 1; frame <= frames; frame++) {
			for (slice = 1; slice <= slices; slice++) {
				Stack.setPosition(channel, slice, frame);;
				applyTfmToSlice(M);
			}
		}
	}
}

function estimateTransform(degree,tblname) {
	// Estimate the alignement for each channel in the image stack
	setBatchMode(true);
	id0 = getImageID();
	Table.create(tblname);
	Stack.getDimensions(width, height, channels, slices, frames);
	if (slices > 1) {
		run("Z Project...", "projection=[Max Intensity]");
	}
	run("Set Measurements...", "  redirect=None decimal=3");
	for (channel = 1; channel <= channels; channel++) {
		selectImage(id0);
		getBeadsLocation(channel);
	}
	setBatchMode(false);
	for (channel = 2; channel <= channels; channel++) {
		Stack.setChannel(channel);
		M = estimateTfm(1,channel,degree);
		cname = "M1" + channel;
		print("Saving model as column " + cname);
		selectWindow(tblname);
		Table.setColumn(cname, M);
		Table.update();
	}
}

function estimateTransformROI(degree,tblname) {
	// Estimate the alignement for each channel in the image stack
	setBatchMode(true);
	id0 = getImageID();
	Table.create(tblname);
	Stack.getDimensions(width, height, channels, slices, frames);
	if (slices > 1) {
		run("Z Project...", "projection=[Max Intensity]");
	}
	getROILocation();
	if (roiManager("count") < 5) {
		exit("Please add more ROIs");
	}
	setBatchMode(false);
	for (channel = 2; channel <= channels; channel++) {
		Stack.setChannel(channel);
		M = estimateTfm(1,channel,degree);
		cname = "M1" + channel;
		print("Saving model as column " + cname);
		selectWindow(tblname);
		Table.setColumn(cname, M);
		Table.update();
	}
}

function estimateTfm(c1,c2,degree) {
	// Estimate the transformation model by loading coordinate from 
	// the result table.
	w = getWidth();
	h = getHeight();
	if (matches(degree,"Affine")) {
		X = loadCoords(c1,3);
	} else {
		X = loadCoords(c1,6);
	}
	//print(mat2str(X, "channel 1: X"));
	if (matcol(X)==0) {
		print("No control points found in channel "+c1); exit();
	}
	
	Y = loadCoords(c2,2);
	if (matcol(X)==0) {
		print("No control points found in channel "+c2); exit();
	}
	//print(mat2str(Y, "channel 2: Y"));

	M = icp(X,Y);

	//print(mat2str(M, "transform"));
	return M;
}

function icp(X,Y) {
	// iterative closest point
	//M = solve(X,Y);
	M = matzeros(X[1],Y[1]);
	matset(M,1,0,1);
	matset(M,2,1,1);	
	for (iter = 0; iter < 10; iter++) {
		MX = matmul(X,M);
		Yi = nnmatch(MX,Y);
		M = solve(X,Yi);
	}
	M = mattranspose(M);
	return M;
}

function nnmatch(X,Y) {
	// for each row of Y[Nx2] find the closest row in X[Nx2]
	Yi = matzeros(X[0],Y[1]);
	for (i = 0; i < X[0]; i++) {
		dstar = 1e6;
		jstar = 0;
		for (j = 0; j < Y[0]; j++) {
			dx = X[2+0+i*2] - Y[2+0+j*2];
			dy = X[2+1+i*2] - Y[2+1+j*2];
			d = dx*dx+dy*dy;
			if (d < dstar) {
				dstar = d;
				jstar = j;
			}
		}
		// copy x and y values in Yi corresponing to the 
		// closest match to Xi
		Yi[2+0+i*2] = Y[2+0+jstar*2];
		Yi[2+1+i*2] = Y[2+1+jstar*2];
	}
	return Yi;
}

function smatch(X,Y) {
	// for each row of Y[Nx2] find the closest row in X[Nx2]
	// use the sinkhorn algorithm to find the best match
	p = mattranspose(X);
	q = mattranspose(Y);
	C = pairwiseDistance(p,q);
	P = sinkhorn_uniform(C,-1);
	mat2img(P,"P");
	Yi = matzeros(X[0],Y[1]);
	for (i = 0; i < X[0]; i++) {
		maxi = -1;
		jstar = 0;
		for (j = 0; j < Y[0]; j++) {
			pij = matget(P,i,j);
			if (pij > maxi) {
				maxi = pij;
				jstar = j;
			}
		}
		// copy x and y values in Yi corresponing to the 
		// closest match to Xi
		Yi[2+0+i*2] = Y[2+0+jstar*2];
		Yi[2+1+i*2] = Y[2+1+jstar*2];
	}
	return Yi;
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
		//print("min c: "+minc+" gamma = "+gamma);
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


function loadCoords(channel, M) {
	// return a NxM matrix loading coordinates from the result table
	N = nResults;
	a = matzeros(N,M);
	k = 0;
	w = getWidth();
	h = getHeight();
	for (i = 0; i < N; i++){
		if (getResult("Channel", i)==channel) {
			x = (getResult("X", i)-w/2)/w*2;
			y = (getResult("Y", i)-h/2)/h*2;
			if (M==6) {
				a[2+M*k] = 1;
				a[2+M*k+1] = x;
				a[2+M*k+2] = y;
				a[2+M*k+3] = x*x;
				a[2+M*k+4] = x*y;
				a[2+M*k+5] = y*y;
			} else if (M==3) {
				a[2+M*k] = 1;
				a[2+M*k+1] = x;
				a[2+M*k+2] = y;
			} else if (M==2) {
				a[2+M*k] = x;
				a[2+M*k+1] = y;
			}
			k = k + 1; 
		}
	}
	a[0] = k-1;
	a = Array.trim(a,2+M*(k-1)+M);
	print("Loaded " + k + " points from channel "+ channel);
	return a;
}

function getBeadsLocation(channel) {
	// Get location of beads in all channels and save the results in a result table
	colors=newArray("red","green","blue","white");
	id = getImageID();
	run("Select None");
	run("Duplicate...", "duplicate channels="+channel);
	run("32-bit");
	run("Gaussian Blur...", "sigma=1");
	id1 = getImageID();
	run("Duplicate...", " ");
	run("Gaussian Blur...", "sigma=2");
	id2 = getImageID();
	imageCalculator("Subtract ", id1, id2);
	getStatistics(area, mean, min, max, std, histogram);
	setThreshold(mean+2*std, max);
	n0 = roiManager("count");
	run("Analyze Particles...", "size=5-Infinity pixel circularity=0.1-1.00 add");
	n1 = roiManager("count");	
	getPixelSize(unit, pixelWidth, pixelHeight);
	selectImage(id);
	for (i = n0; i < n1; i++){		
		roiManager("select",i);
		setResult("Channel",i, channel);
		xm = getValue("XM")/pixelWidth;
		ym = getValue("YM")/pixelWidth;
		Stack.setChannel(channel);
		pos = fitBeadXY(xm,ym,2);		
		if (!(isNaN(pos[0]) || isNaN(pos[0]))) {
			xm = pos[0];
			ym = pos[1];
		}
		setResult("X",i,xm);
		setResult("Y",i,ym);		
		Overlay.addSelection(colors[(channel-1)%colors.length]);
		Overlay.show;
		updateResults();
	}
	updateResults();
	selectImage(id1);close();
	selectImage(id2);close();
}

function fitBeadXY(x0,y0,s) {
	bg = 1e9;
	for (dy = -s; dy <= s; dy+=0.5) {
		for (dx = -s; dx <= s; dx+=0.5) {
			bg = minOf(bg, getPixel(x0+dx, y0+dy));
		}
	}
	for (i = 0; i < 10; i++) {			
		sv = 0;
		x1 = 0;
		y1 = 0;
		for (dy = -s; dy <= s; dy+=0.5) {
			for (dx = -s; dx <= s; dx+=0.5) {
				d = (dx*dx+dy*dy)/(9*s*s);			
				v = maxOf(0, getPixel(x0+dx, y0+dy)-bg) * exp(-0.5*d);				
				x1 += (x0+dx) * v;
				y1 += (y0+dy) * v;
				sv += v;
			}
		}		
		x0 = x1/sv;
		y0 = y1/sv;
		wait(1);		
	}	
	return newArray(x0,y0);
}

function getROILocation() {
	// Get location of rois in all channels and save the results in a result table
	run("Select None");
	n = roiManager("count");
	getPixelSize(unit, pixelWidth, pixelHeight);
	for (i = 0;i < n; i++){		
		roiManager("select",i);
		Stack.getPosition(channel, slice, frame);
		setResult("Channel",i, channel);
		setResult("X",i,getValue("XM")/pixelWidth);
		setResult("Y",i,getValue("YM")/pixelHeight);
		Overlay.addSelection();
		Overlay.show;
		updateResults();
	}
	updateResults();	
}


function applyTfm(filename) {
	// apply a tranform from a file to the all dataset
	Stack.getDimensions(width, height, channels, slices, frames);
	for (c = 1; c < channels; c++) {
		tfm = loadTfm(filename, c);
		for (t = 1; t < frames; t++) {
			for (z = 1; z < slices; z++) {
				Stack.setPosition(channel, slice, frame);
				applyTfmToSlice(tfm);
			}
	}	
}

function applyTfmToSlice(M) {
	// apply a transform M to the current image slice
	w = getWidth;
	h = getHeight;
	buf = newArray(w*h);
	if (M[1] == 6) {
		for (y = 0; y < h; y++) {
			for (x = 0; x < w; x++) {
				xi = (x - w/2)/w*2;
				yi = (y - h/2)/h*2;
				xt = M[2] + M[3] * xi + M[4] * yi +  M[5] * xi*xi + M[6] * xi*yi + M[7]*yi*yi;
				yt = M[8] + M[9] * xi + M[10] * yi + M[11] * xi*xi + M[12] * xi*yi + M[13]*yi*yi;
				xt = xt*w/2 + w/2;
				yt = yt*h/2 + h/2;
				if (xt >= 0 && yt >= 0 && xt < w && yt < h) {
					i = x + getWidth * y;
					buf[i] = getPixel(xt,yt);
				}
			}
		}
	} else {
		for (y = 0; y < h; y++) {
			for (x = 0; x < w; x++) {
				xi = (x - w/2)/w*2;
				yi = (y - h/2)/h*2;
				xt = M[2] + M[3] * xi + M[4] * yi;
				yt = M[5] + M[6] * xi + M[7] * yi;
				xt = xt*w/2 + w/2;
				yt = yt*h/2 + h/2;
				if (xt >= 0 && yt >= 0 && xt < w && yt < h) {
					i = x + getWidth * y;
					buf[i] = getPixel(xt,yt);
				}
			}
		}
	}
	for (y = 0; y < h; y++) {
		for (x = 0; x < w; x++) {
			i = x + getWidth * y;
			setPixel(x,y,buf[i]);
		}
	}
}

function image2Array(w,h) {
	// Conver the image to an 1D array
	A = newArray(w*h);
	for (y = 0; y < h; y++) {
		for (x = 0; x < w; x++) {
			i = x + getWidth * y;
			A[i] = getPixel(x, y);
		}
	}
	return A;
}



/////////////////////////////////////////////////////////
// Basic Linear Algebra
////////////////////////////////////////////////////////
function assert(cond, msg) {
	// check that condition is true
	if (!(cond)) {
		exit(msg);
	}
}

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
		for (m = 0; m < M; m++) {
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
	//print(mat2str(A,"A"));
	//print(mat2str(B,"B"));
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


function testmat() {
	// a serie of test for matrices
	A = matnew(2,2,newArray(1,0,0,1));
	B = matnew(2,2,newArray(2,1,0,1));
	print(mat2str(A,"A"));
	print(mat2str(B,"B"));
	C = matadd(A,B,1);
	print(mat2str(C,"A + B"));
	C = matmul(A,B);
	print(mat2str(C,"A x B"));
	C = mattranspose(B);
	print(mat2str(C,"Bt"));
	C = solve(A,B);
	print(mat2str(C,"A\\B"));
	C = matrep(A,2,3);
	print(mat2str(C,"repmat(A,2,3)"));
	C = matsum(A,0);
	print(mat2str(C,"matsum(A,0)"));
	C = matmean(A,0);
	print(mat2str(C,"matmean(A,0)"));
}

/// END Of Basic Linear Algebra routines
