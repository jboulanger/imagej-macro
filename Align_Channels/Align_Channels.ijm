#@String(choices={"Estimate Beads","Estimate ROI","Apply","Show beads","Test"}) mode
#@String(choices={"Affine","Quadratic"}) model
#@String(label="table name",value="Models.csv") tblname

/*
 * Estimate and apply channel alignement from a bead calibration image
 *
 * Channels are aligned to the 1st channel
 *
 * Estimate transform:
 *  - Open an image with beads and save the model as a csv file
 *  - If the image is a 3d stack a mip will be performed first
 *
 * Apply transform:
 *   - Parameters are loaded from a csv table
 *   - All planes of the hyperstacks are aligned
 *
 * If there is no opened image a test image is generated.
 *
 *  Toubleshouting
 *
 *  	1. The macro will work on small images (tested on 1024x1024) and is relatively slow.
 *  	2. Some error may occur, try closing the results and models.csv tables.
 *
 *
 *  TODO: improve robustness to outliers, use multi-frame data
 *
 * Jerome Boulanger 2021-2022 for Sina, Anna & Ross
 */

run("Set Measurements...", "  redirect=None decimal=12");

if (nImages==0) {
	print("No image opened, generating a test image with a quadratic model.");
	generateTestImage();
	model="Quadratic";
	estimateTransform(model,tblname);
	run("Select None");
	run("Duplicate...","title=[test aligned] duplicate");
	Overlay.remove();
	applyTransform(tblname);
	measureTestError(getImageID());
	exit();
}

if (matches(mode, "Estimate Beads")) {
	print("Estimate transfrom from beads");
	estimateTransform(model,tblname);
} else if (matches(mode, "Estimate ROI")) {
	print("Estimate transfrom from ROIs");
	estimateTransformROI(model,tblname);
} else if (matches(mode,"Show beads")) {
	showBeads();
} else if (matches(mode,"Apply")){
	print("Apply transform");
	applyTransform(tblname);
} else {
	print("Test");
	run("Close All");
	generateTestImage();
	model="Quadratic";
	estimateTransform(model,tblname);
	run("Select None");
	run("Duplicate...","title=[test aligned] duplicate");
	Overlay.remove();
	applyTransform(tblname);
	measureTestError(getImageID());
}
print("Done");

function generateTestImage() {
	w = 500;
	h = 500;
	newImage("test", "8-bit composite-mode", w, h, 2, 1, 1);
	N = 10;
	a = 2e-2;
	b = 1e-3;
	c = 1e-4;
	M = newArray(a*random,1+b*random,b*random,c*random,c*random,c*random,   a*random,b*random,1+b*random,c*random,c*random,c*random);
	//Array.print(M);
	rho = 0.9;
	for (k = 0; k < N; k++) {
		x = w/8 + 3/4*w * random;
		y = h/8 + 3/4*h * random;
		Stack.setChannel(1);
		if (random < rho) {
			drawGaussian(x,y);
		}
		xi = 2 * x / w - 1;
		yi = 2 * y / h - 1;
		xt = M[0] + M[1] * xi + M[2] * yi +  M[3] * xi*xi + M[4] * xi*yi + M[5]*yi*yi;
		yt = M[6] + M[7] * xi + M[8] * yi + M[9] * xi*xi + M[10] * xi*yi + M[11]*yi*yi;
		xt = w * (xt + 1) / 2;
		yt = h * (yt + 1) / 2;
		Stack.setChannel(2);
		if (random < rho) {
			drawGaussian(xt,yt);
		}
	}
}

function drawGaussian(xc,yc) {
	x0 = round(xc);
	y0 = round(yc);
	for (y = y0-3; y <= y0+3; y++) {
		for (x = x0-3; x <= x0+3; x++) {
			dx = x - xc;
			dy = y - yc;
			I = getPixel(x, y);
			setPixel(x, y, I + 128 * exp(-0.25*(dx*dx+dy*dy)));
		}
	}
}

function measureTestError(id0) {
	selectImage(id0);
	run("Select None");
	Stack.setChannel(1);
	run("Duplicate...", "title=c1 channels=1");
	run("32-bit");
	maxi = getValue("Max");
	id1 = getImageID();
	selectImage(id0);
	Stack.setChannel(2);
	run("Duplicate...", "title=c2 channels=2");
	run("32-bit");
	id2 = getImageID();
	imageCalculator("Subtract", id1, id2);
	run("Square");
	run("Select None");
	rmse = sqrt(getValue("Mean"));
	print("RMSE : " + rmse);
	print("PSNR  : " + 20 * log(maxi / rmse) +"dB");
	selectImage(id1); close();
	selectImage(id2); close();
	return rmse;
}

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
	Overlay.remove;
	setBatchMode(true);
	if (isOpen("Results")) {selectWindow("Results");run("Close");}
	id0 = getImageID();
	Table.create(tblname);
	Stack.getDimensions(width, height, channels, slices, frames);
	if (slices > 1) {		
		print("Max intensity projection");
		run("Z Project...", "projection=[Max Intensity]");
	}

	run("Set Measurements...", "  redirect=None decimal=6");
	for (frame = 1; frame <= frames; frame++) {
		for (channel = 1; channel <= channels; channel++) {
			selectImage(id0);
			getBeadsLocation(channel, frame);
		}
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
	Overlay.remove;
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

function icpold(X,Y) {
	M = matzeros(X[1],Y[1]);
	matset(M,1,0,1);
	matset(M,2,1,1);
	print(mat2str(M,"M"));
	
	for (iter = 0; iter < 2; iter++) {
		MX = matmul(X,M);
		print(mat2str(MX,"MX"));
		I = nnmatch(MX,Y);
		XI = matpermuterow(X, matextractcolumn(I, 0));
		YI = matpermuterow(Y, matextractcolumn(I, 1));
		print(mat2str(XI,"XI"));
		print(mat2str(YI,"YI"));
		M = solve_from(XI,YI,M,100);		
		print(mat2str(M,"M"));
	}
	
	setColor("yellow");	drawMatch(XI, YI);
	M = mattranspose(M);
	return M;
}

function icp(X,Y) {
	// iterative closest point
	// X : coordinate matrix
	// Y : coordinate matrix
	// return the transformation matrix
	n = X[0];
	M = matzeros(X[1],Y[1]);
	matset(M,1,0,1);
	matset(M,2,1,1);
	
	e1 = 1e6;

	for (iter = 0; iter < 2; iter++) {

		MX = matmul(X,M);		
		I = smatch(MX,Y);		
		XI = matpermuterow(X, matextractcolumn(I, 0));
		YI = matpermuterow(Y, matextractcolumn(I, 1));		
		M = solve_from(XI,YI,M,10000,1e-9);
		MXI = matmul(XI,M);		
		E = errorDistance(MXI,YI);
		e0 = e1;
		e1 = matmean(E,-1);
		
		/*
		for (inner = 0; inner < 1; inner++) {
			E = errorDistance(MXI,YI);						
			J = matbelowthreshold(E, 1.5*e1);
			XIs = matsubsetrow(XI, J);
			YIs = matsubsetrow(YI, J);			
			print("number of points " + XIs.length);
			M = solve(XIs, YIs);
			MXI = matmul(XI,M);
		}
		
		M = solve(XI, YI);
		MXI = matmul(XI,M);
		E = errorDistance(MXI,YI);
		e0 = e1;
		e1 = matmean(E,-1);				
		print("error : " + e1 + ", relative variation: "+ (e1-e0) / (e1+e0));
		if (iter > 0 && abs(e1-e0)/(e1+e0)<1e-12) { break++; }
		*/
	}
	//MXs = matmul(Xs,M);
	M = mattranspose(M);
	setColor("yellow");	drawMatch(XI, YI);
	return M;
}

function drawMatch(X,Y) {
	// draw a line between each rows of X and Y
	getDimensions(w, h, channels, slices, frames);
	for (i = 0; i < matrow(X); i++) {		
		if (matcol(X)==2) {
			makeLine(w*(matget(X,i,0)+1)/2, h*(matget(X,i,1)+1)/2, w*(matget(Y,i,0)+1)/2, h*(matget(Y,i,1)+1)/2);
		} else {
			makeLine(w*(matget(X,i,1)+1)/2, h*(matget(X,i,2)+1)/2, w*(matget(Y,i,0)+1)/2, h*(matget(Y,i,1)+1)/2);
		}
		Roi.setStrokeWidth(0.1);
		Overlay.addSelection("yellow");		
		//Overlay.show();
	}
	run("Select None");
}

function nnmatch(X,Y) {
	// for each row of Y[Nx2] find the closest row in X[Nx2]
	n = matrow(X);
	m = matrow(Y);
	I = matzeros(n*m, 2);
	k = 0;	
	for (i = 0; i < X[0]; i++) {
		dstar = 1e6;
		jstar = 0;
		for (j = 0; j < Y[0]; j++) {
			dx = X[2+0+i*2] - Y[2+0+j*2];
			dy = X[2+1+i*2] - Y[2+1+j*2];
			d = sqrt(dx*dx+dy*dy);
			if (d < dstar) {
				dstar = d;
				jstar = j;
			}
		}
		// copy x and y values in Yi corresponing to the
		// closest match to Xi
		if (dstar < 0.05) {
			matset(I,k,0,i);
			matset(I,k,1,jstar);
			k = k + 1;				
		}		
	}
	I[0] = k;
	return I;
}

function smatch(X,Y) {
	/* Sinkhorn match
	 * X : [N,2] coordinates
	 * Y : [N,2] coordinates
	 */
	n = matrow(X);
	m = matrow(Y);
	C = computeAssignmentCostMatrix(X, Y);
	//C = pairwiseDistance(X,Y);
	//id=getImageID();mat2img(C,"C");selectImage(id);
	P = sinkhorn_uniform(C,0.05);
	pmax = matmax(P,-1);
	//id = getImageID(); mat2img(P,"P"); wait(100); selectImage(id);

	// find candidate matches
	I = matzeros(n*m, 2);	
	k = 0;	
	for (i = 0; i < n; i++) {						
		// find max along j
		jmaxi = -1;
		maxi = -1;
		for (j = 0; j < m; j++) {
			val = matget(P,i,j);
			if (val > maxi && val > pmax/2) {
				jmaxi = j;
				maxi = val;
			}
		}
		if (jmaxi == -1) {continue;}
		
		// check if the also max along i
		imaxi = -1;
		maxi = -1;
		for (ii = 0; ii < n; ii++){
			val = matget(P,ii,jmaxi);
			if (val >= maxi) {
				imaxi = ii;
				maxi = val;
			}
		}
				
		// record the pair
		if (imaxi == i && jmaxi >= 0 ) {
			matset(I,k,0,i);
			matset(I,k,1,jmaxi);
			k = k + 1;				
		}
	}	
	// set the number of row to k
	I[0] = k;	
	return I;
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
		minc = matmin(C,-1);
		gamma = maxOf(1e-5, 200*minc);		
		print("\n*** min c: " + minc + " gamma = " + gamma + " ***\n");
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
	for (iter = 0; iter < 50; iter++) {
		xib = matmul(xi,b);
		a = mat_div(p,xib);
		xita = matmul(xit,a);
		b = mat_div(q,xita);
	}
	// P = diag(a)*xi*diag(b)
	P = matzeros(n,m);
	for (i = 0; i < n; i++) {
		for (j = 0; j < m; j++) {
			P[2+j+m*i] = xi[2+j+m*i] * a[2+i] * b[2+j];
		}
	}
	/*
	a = matdiag(a);
	b = matdiag(b);
	P = matmul(xi, b);
	P = matmul(a, P);
	*/
	return P;
}

function meanErrorDistance(X,Y) {
	// Compute the error distance between matching X and Y
	// X : [n,2] coordinate matrix
	// Y : [m,2] coordinate matrix
	// return the average distance
	//print("X ["+X[0]+","+X[1]+"]");
	//print("Y ["+Y[0]+","+Y[1]+"]");
	n = matrow(X);
	D = 2;
	e = 0;
	for (i = 0; i < n; i++) {
		s = 0;
		for (k = 0; k < D; k++) {
			d = matget(X,i,k) - matget(Y,i,k);
			s += d*d;
		}
		e += sqrt(s);
	}
	e = e/n;
	return e;
}

function errorDistance(X,Y) {
	// Compute the error distance between matching X and Y
	// X : [n,2] coordinate matrix
	// Y : [m,2] coordinate matrix
	// return error
	//print("X ["+X[0]+","+X[1]+"]");
	//print("Y ["+Y[0]+","+Y[1]+"]");
	n = matrow(X);
	D = 2;
	E = matzeros(n,1);
	for (i = 0; i < n; i++) {
		s = 0;
		for (k = 0; k < D; k++) {
			d = matget(X,i,k) - matget(Y,i,k);
			s += d*d;
		}
		matset(E,i,0,sqrt(s));
	}
	return E;
}

function pairwiseDistance(p,q) {
	// Compute the pairwise distance between p and q
	// P : coordinate matrix [n,2]
	// Q : coordinate matrix [m,2]
	// return distance [n x m] matrix
	//print("p ["+p[0]+","+p[1]+"]");
	//print("q ["+q[0]+","+q[1]+"]");
	n = matrow(p);
	m = matrow(q);
	D = 2;
	a = matzeros(n,m);
	for (i = 0; i < n; i++) {
		for (j = 0; j < m; j++) {
			s = 0;
			for (k = 0; k < D; k++) {
				d = matget(p,i,k) - matget(q,j,k);
				s += d*d;
			}
			matset(a,i,j,s);
		}
	}
	return a;
}

function computeAssignmentCostMatrix(p,q) {
	/* Cost matrix with missing points
	 *
	 * P : coordinate matrix [n,2]
	 * Q : coordinate matrix [m,2]
	 * return distance [n+m x m+n] matrix
	 */
	n = matrow(p);
	m = matrow(q);
	D = 2;
	a = matzeros(n+m,m+n);	
	cmax = 0;
	for (i = 0; i < n; i++) {
		for (j = 0; j < m; j++) {
			s = 0;
			for (k = 0; k < D; k++) {
				d = matget(p,i,k) - matget(q,j,k);
				s += d*d;
			}
			matset(a,i,j, sqrt(s));
			cmax += sqrt(s);
		}
	}
	cmax = matmax(a,-1);
	//cmax = 0.2;//matmean(a,-1);	
	for (i = 0; i < n; i++) {
		for (j = 0; j < m; j++) {
			if (matget(a,i,j) >	cmax) {
				matset(a,i,j, cmax);
			}
		}
	}
	
	for (i = 0; i < n; i++) {
		for (j = m; j < m+n; j++) {
			matset(a,i,j,cmax);
		}
	}
	for (i = n; i < n+m; i++) {
		for (j = 0; j < m; j++) {
			matset(a,i,j,cmax);
		}
	}
	for (i = n; i < n+m; i++) {
		for (j = m; j < n+m; j++) {
			matset(a,i,j,2*cmax);
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
	a[0] = k;
	a = Array.trim(a,2+M*(k)+M);
	print("Loaded " + k + " points from channel "+ channel);
	return a;
}

function getBeadsLocation(channel, frame) {
	// Get location of beads in all channels and save the results in a result table
	
	colors = newArray("red","green","blue","white");
	id = getImageID();
	run("Select None");
	run("Duplicate...", "duplicate channels="+channel+" frames="+frame);
	run("32-bit");
	run("Gaussian Blur...", "sigma=2");
	id1 = getImageID();
	run("Duplicate...", " ");
	run("Gaussian Blur...", "sigma=6");
	id2 = getImageID();
	imageCalculator("Subtract ", id1, id2);
	getStatistics(area, mean, min, max, std, histogram);
	setThreshold(mean+4*std, max);
	n0 = roiManager("count");
	run("Analyze Particles...", "size=10-Infinity pixel circularity=0.1-1.00 add");
	n1 = roiManager("count");
	getPixelSize(unit, pixelWidth, pixelHeight);
	selectImage(id);
	for (i = n0; i < n1; i++){
		roiManager("select",i);
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
		setResult("Channel",i, channel);
		setResult("Frame",i, frame);
		//Overlay.drawEllipse(xm-1, ym-1, 2, 2);
		//Overlay.drawEllipse(100,100,2,2);
		//Overlay.show;
		makeOval(xm-0.5, ym-0.5, 2, 2);
		Roi.setStrokeColor(colors[(channel-1)%colors.length]);
		Roi.setStrokeWidth(0.1);
		Overlay.addSelection();
		
		//Overlay.setStrokeColor(colors[(channel-1)%colors.length]);
		//Overlay.addSelection(colors[(channel-1)%colors.length]);
		//Overlay.show;
		updateResults();
	}
	updateResults();
	selectImage(id1);close();
	selectImage(id2);close();
	run("Select None");
	
}

function fitBeadXY(x0,y0,s) {
	bg = 1e9;
	for (dy = -s; dy <= s; dy+=0.5) {
		for (dx = -s; dx <= s; dx+=0.5) {
			bg = minOf(bg, getPixel(x0+dx, y0+dy));
		}
	}
	for (i = 0; i < 30; i++) {
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
	for (i = 0; i < n; i++){
		roiManager("select",i);
		Stack.getPosition(channel, slice, frame);
		setResult("X",i,getValue("XM")/pixelWidth);
		setResult("Y",i,getValue("YM")/pixelHeight);
		setResult("Channel",i, channel);
		setResult("Frame",i, frame);
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

function matcheck(A, name) {
	// check if the size of the matrix are consistent
	if (A.length != 2+A[0]*A[1]) {
		print("Matrix "+name+" is not consistent");
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
	// return a reshaped matrix
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
	//assert((A.length == 2+N*M),"mattranspose: size of the array does not match matrix dimensions");

	B = matzeros(M,N);
	for (n = 0; n < N; n++) {
		for (m = 0; m < M; m++) {
			B[2+n+m*N] = A[2+m+n*M];
		}
	}
	return B;
}

function matshape2str(A,name) {
	if (A.length == 0){
		str = name + "is an empty matrix";
	} else {
		str = name + " ["+A[0]+"x"+A[1]+"] size:"+A.length+"\n";
	}
	return str;
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

function matmin(A,dim) {
	// average along one dimension
	N = A[0];
	M = A[1];
	if(dim==0) {
		S = matzeros(1,M);
		for (n = 0; n < N; n++) {
			for (m = 0; m < M; m++) {
				S[2+m] = minOf(S[2+m],A[2+m+n*M]);
			}
		}
	} else if (dim==1){
		S = matzeros(N,1);
		for (n = 0; n < N; n++) {
			for (m = 0; m < M; m++) {
				S[2+n] = minOf(S[2+n],A[2+m+n*M]);
			}
		}
	} else {
		S = A[2];
		for (n = 0; n < N; n++) {
			for (m = 0; m < M; m++) {
				S = minOf(S,A[2+m+n*M]);
			}
		}
	}
	return S;
}

function solve(A,B) {
	// Solve A*X=B by iterating X=X-At(AX - B)
	// A[NxM] x X[MxK] = B [NxK]
	
	N = A[0];
	M = A[1];
	K = B[1];
	At = mattranspose(A);
	AtA = matmul(At,A);
	AtB = matmul(At,B);
	X = matzeros(M,K);
	a = matnorm(AtA);
	for (iter = 0; iter < 10000; iter++) {
		AtAX = matmul(AtA,X);
		dX = matadd(AtAX,AtB,-1);
		X = matadd(X, dX, -1/(2*a));
		if (matmax(dX,-1) < 1e-9) {break;}
	}
	return X;
}

function solve_from(A,B,X0,niter,eps) {
	// Solve A*X=B by iterating X=X-At(AX - B)
	// A[NxM] x X[MxK] = B [NxK]
	
	N = A[0];
	M = A[1];
	K = B[1];
	At = mattranspose(A);
	AtA = matmul(At,A);
	AtB = matmul(At,B);
	X = X0;
	a = matnorm(AtA);
	for (iter = 0; iter < niter; iter++) {
		AtAX = matmul(AtA,X);
		dX = matadd(AtAX,AtB,-1);
		X = matadd(X, dX, -1/(2*a));
		delta = maxOf(abs(matmin(dX,-1)), abs(matmax(dX,-1)));
		if (delta < eps) {break;}
	}
	print("Stopping after " + iter + "/" + niter+ ", delta: " + delta);
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

function matabovethreshold(A,t) {
	/* Threshold matrix A
	 *  returns a matrix with 1 if Aij> t
	 */
	B = matzeros(matrow(A),matcol(A));
	for (i=2;i<A.length;i++) {
		if (A[i] > t) {B[i]=1;}
	}
	return B;
}

function matbelowthreshold(A,t) {
	/* Threshold matrix A
	 *  returns a matrix with 1 if Aij> t
	 */
	B = matzeros(matrow(A),matcol(A));
	for (i=2;i<A.length;i++) {
		if (A[i] < t) {B[i]=1;}
	}
	return B;
}

function matsubsetrow(A,I) {
	// Extract rows from A where I is 1
	n = matrow(A);
	m = matcol(A);
	if (I[0] != A[0]) {
		print("error in matsubsetrow, I has "+I[0]+" rows and A has "+A[0]" rows");
	}
	k = 0;
	B = Array.copy(A);
	for (i = 0; i < n; i++) {
		if (I[i] > 0) {
			for (j = 0; j < m; j++) {
				B[2+j+m*k] = A[2+j+m*i];
				//matset(B,k,j,matget(A,i,j));
			}
			k = k+1;
		}
	}
	B = Array.trim(B, 2+m*(k));
	B[0] = k;
	return B;
}

function matpermuterow(A,I) {
	n = matrow(I);
	m = matcol(A);
	B = matzeros(n, matcol(A));
	for (i = 0; i < n; i++) {
		for (j = 0; j < m; j++) {
			matset(B, i, j, matget(A, matget(I,i,0), j));
		}
	}
	return B;
}

function matextractcolumn(A, j) {
	/* Extract the jth column from A*/
	n = matrow(I);
	B = matzeros(n, 1);
	for (i = 0; i < n; i++) {
		matset(B, i, 0, matget(A, i, j));
	}
	return B;
}

function matsampler(n,m,rho) {
	/* Generate a random sampling matrix
	 * n: number of rows
	 * m : number of colmns
	 * rho: proba [0-1] of being 1
	 */
	I = matzeros(n,m);
	for (i=2;i<n*m;i++) {
		if (random < rho) {
			I[i] = 1;
		}
	}
	return I;
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
