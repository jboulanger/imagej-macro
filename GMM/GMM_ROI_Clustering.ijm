#@String(label="Measurements", value="AR,Feret") featstr
#@String(label="Classes", value="Class1,Class2") classstr
/*
 * GMM ROI Clustering
 * 
 * Clustering of ROI based on features using a Gaussian Mixture Model
 * 
 * Usage:
 * - Add ROI in the ROI manager (Analyse particle for example)
 * - Run the macro and select measurements and classes
 * 
 * To do:
 * - multi-start initialization
 * 
 */
features = split(featstr, ",");
classes = split(classstr, ",");

if (isOpen("ROI Manager")) {
	id = getImageID();
	if(isOpen("plot")) { selectWindow("plot"); run("Close"); }
	if(isOpen("gmm")) { selectWindow("gmm"); run("Close"); }	
	npts = roiManager("count");
	K = gmmClusterROI(features, classes);
	selectImage(id);
	X = getMeasureValueForAllROI(features);
	stats = gmmComputeStatsByClass(X,K,npts,features.length,classes.length);
	gmmPrintStats(stats, features, classes);
	gmmSaveStatsToTable("Analysis.csv", "test.tif","test",stats,features,classes);
} else {
	print("Demo mode");
	run("Close All");
	npts = 1000;
	nclass = classes.length;
	nfeat = features.length;
	C0 = createClasses(npts,nclass);
	X = createData(C0,npts,nfeat,nclass);
	plotData(X,C0,npts,features,nclass);
	start = getTime();
	C = gmmFit(X, npts, nfeat, nclass, 0.7, false,30);
	print("Elapsed time "+ (getTime()-start)/1000+" s");
	plotData(X,C,npts,features,nclass);
	computeClassifierPerf(C0,C);
	res = gmmComputeStatsByClass(X,C,npts,nfeat,nclass);
	//gmmPrintStats(stats, features, newArray("A","B"));
}

function createClasses(n,nclass) {
	// create classes
	C = newArray(n);
	for (k = 0; k < nclass; k++) {
		p = 1 / nclass;
		for (j = 0; j < n; j++) {
			if (random < p) { 
				C[j] = k;
			}
		}
	}
	return C;
}

function createData(K,npts,nfeat,nclass) {
	// create the synthetic test data
	X = newArray(npts*nfeat);
	C = newArray(nfeat*nclass);
	S = newArray(nfeat*nclass);
	for (k = 0; k < nclass; k++) {
		for (i = 0; i < nfeat; i++) {
			if (k > 0) {
				C[i+nfeat*k] = 7 * random;
			}
			S[i+nfeat*k] = 0.5 + 1 * random;	
		}
	}
	for (j = 0; j < npts; j++) {
		k = K[j];
		for (i = 0; i < nfeat; i++) {
			X[i+nfeat*j] = C[i+nfeat*k] + S[i+nfeat*k] * random("gaussian");
		}
	}
	return X;
}

function computeClassifierPerf(C0,C) {
	fp = 0;
	fn = 0;
	tp = 0;
	tn = 0;
	for (j=0;j<npts;j++) {
		if (C0[j]==0 && C[j]==1) {fp++;}
		if (C0[j]==1 && C[j]==0) {fn++;}
		if (C0[j]==0 && C[j]==0) {tn++;}
		if (C0[j]==1 && C[j]==1) {tp++;}
	} 
	print("=== Performances ===");
	print("False positive rate "+ fp / (fp+tn) );
	print("False negative rate "+ fn / (tp+fn) );
	print("Specificity (1-FPR) "+ 1 - fp / (fp+tn) );
	print("Sensitivity (1-FNR) "+ 1 - fn / (tp+fn) );
	print("Precision (TP/Detected)"+ tp / (tp+fp) );
}



// start of clustering code

function gmmClusterROI(features, classes) {
	/* Cluster ROI uing features
	 * Parameter
	 * features : array of features (string)
	 * Result
	 * C array with the index the class of each ROI
	 * Example
	 * gmmClusterROI(newArray("AR","Area"));
	 */
	X = getMeasureValueForAllROI(features);
	for (i=0;i<X.length;i++) {X[i] = log(X[i]);}
	npts = roiManager("count");
	nfeat = features.length;
	nclass = classes.length;
	//K = newArray(npts);
	K = gmmFit(X, npts, nfeat, nclass, 1.0, true, 30);
	colors = newArray("red","blue","green","yellow");
	for (j = 0; j < npts; j++) {
		roiManager("select", j);
		Roi.setStrokeColor(colors[K[j]]);
	}
	plotData(X,K,npts,features,nclass);
	res = gmmComputeStatsByClass(X,K,npts,nfeat,nclass);
	return K;
}

function plotData(X,C,npts,features,nclass) {
	// plot the 2 first components of the data by classes
	nfeat = features.length;
	Plot.create("plot", features[0], features[1]);
	colors = newArray("red","blue","green","yellow");
	markers = newArray("circle","circle","circle");
	xmin = 1e3; xmax = -1e3; ymin = 1e3; ymax = -1e3;
	for (k = 0; k < nclass; k++) {
		xvalues = newArray(npts);
		yvalues = newArray(npts);
		l = 0;
		for (j = 0; j < npts; j++) {
			if (C[j]==k) {
				xvalues[l] = X[nfeat*j];
				yvalues[l] = X[1+nfeat*j];
				xmin = minOf(xmin, xvalues[l]);
				xmax = maxOf(xmax, xvalues[l]);
				ymin = minOf(ymin, yvalues[l]);
				ymax = maxOf(ymax, yvalues[l]);
				l = l + 1;
			}
		}
		xvalues = Array.trim(xvalues,l);
		yvalues = Array.trim(yvalues,l);
		Plot.setColor(colors[k]);
		Plot.add(markers[k], xvalues, yvalues);
	}
	Plot.setLimits(xmin,xmax,ymin,ymax);
	Plot.show();
}

function getMeasureValueForAllROI(features) {
	// return an array (n,m) containing the values of measurements defined by the
	// array measures for each ROI of the ROI manager.
	npts = roiManager("count");
	nfeat = features.length;
	x = newArray(npts * nfeat);
	for (j = 0; j < npts; j++) {
		roiManager("select", j);
		 List.setMeasurements();
		for (i = 0; i < nfeat; i++) {
			x[i+nfeat*j] = List.getValue(features[i]);
		}
	}
	return x;
}

function gmmFit(X, npts, nfeat, nclass, bsize, display, niter) {
	// gmm clustering 
	P = newArray(nclass * npts); // probablity
	C = newArray(nclass * nfeat); // centers
	S = newArray(nclass * nfeat); // std
	// initialization
	gmmInitialize(X,C,S,P,npts,nfeat,nclass);
	if (display) {
		newImage("gmm", "32-bit composite-mode", 100, 100, nclass, 1, 1);
		visu = getImageID();
	} else {
		visu = 1;
	}
	for (iter = 0; iter < niter; iter++) {
		gmmUpdateProbabilities(X,C,S,P,npts,nfeat,nclass,1);
		gmmUpdateCentersMean(X,C,S,P,npts,nfeat,nclass,bsize);
		gmmUpdateImage(visu,C,S,nclass);
	}
	gmmUpdateProbabilities(X,C,S,P,npts,nfeat,nclass,1);
	K = gmmProba2class(P,npts,nclass);
	return K;
}


function gmmInitialize(X,C,S,P,npts,nfeat,nclass) {
	// initialize 2 classes with large/small values for each features
	for (j = 0; j < npts; j++) {
		if (X[1+nfeat*j] > 2  - 1 * X[nfeat*j]) {
			P[nclass*j] = 0;
			P[1+nclass*j] = 1;
		} else {
			P[nclass*j] = 1;
			P[1+nclass*j] = 0;		
		}
	}
	gmmUpdateCentersMean(X,C,S,P,npts,nfeat,nclass,1);
}

function gmmProba2class(P,npts,nclass) {
	// argmax P along class axis
	C = newArray(npts);
	for (j = 0; j < npts; j++) {
		pmax = 0;
		kmax = 0;
		for (k = 0; k<nclass;k++) {
			if (P[k+j*nclass] > pmax) {
				pmax = P[k+j*nclass];
				kmax = k;
			}
		}
		C[j] = kmax;
	}
	return C;
}

function gmmUpdateProbabilities(X,C,S,P,npts,nfeat,nclass,bsize) {
	/** Update the probabilies for each point 
	 *  Parameter
	 *  	X : features (nfeat, npts)
	 *  	C : centers (nfeat, nclass)
	 *  	S : standard deviation (nfeat, nclass)
	 *  	P : probabilities (nclass,npts)
	 *  	npts : number pf data points
	 *  	nfeat : number of features
	 *  	nclass : number of classes
	 *  Return
	 *  	updates C and S inplace
	 */

	// initialize to 0
	for (i = 0; i < P.length; i++) {
		P[i] = 0;
	}
	// compute for each class and point the squared distance of the point Xj to the center Ck
	for (k = 0; k < nclass; k++) {
		for (j = 0; j < npts; j++)  if (random < bsize){
			for (i = 0; i < nfeat; i++) {
				delta = (X[i+nfeat*j] - C[i+nfeat*k]) / S[i+nfeat*k];
				P[k+nclass*j] = P[k+nclass*j] + delta*delta; 
			}
		}
	}
	// take the Gaussian kernel on the distance
	for (i = 0; i < P.length; i++) {
		P[i] = exp(-0.5*P[i]);
	}
	// Ensure sum P is 1 for each point 
	for (j = 0; j < npts; j++) {
		s = 0;
		for (k = 0; k < nclass; k++) {
			s += P[k+nclass*j];
		}
		for (k = 0; k < nclass; k++) {
			P[k+nclass*j] = P[k+nclass*j] / s;
		}
	}
}

function gmmUpdateCentersMean(X,C,S,P,npts,nfeat,nclass,bsize) {
	/** Update the center and standard deviation of each cluster 
	 *  Parameter
	 *  	X : features (nfeat, npts)
	 *  	C : centers (nfeat, nclass)
	 *  	S : standard deviation (nfeat, nclass)
	 *  	P : probabilities (nclass,npts)
	 *  	npts : number pf data points
	 *  	nfeat : number of features
	 *  	nclass : number of classes
	 *  Return
	 *  	updates C and S inplace
	 */
	// Initialize the center accumulator to 0
	for (i = 0; i < C.length; i++) {
		C[i] = 0;
		S[i] = 0;
	}
	// and the partition (sum of P)
	Z = newArray(nclass * nfeat);
	for (j = 0; j < npts; j++) if (random < bsize) {
		for (i = 0; i < nfeat; i++) {
			for (k = 0; k < nclass; k++) {
				C[i+nfeat*k] = C[i+nfeat*k] + P[k+nclass*j] * X[i+nfeat*j];
				S[i+nfeat*k] = S[i+nfeat*k] + P[k+nclass*j] * X[i+nfeat*j] * X[i+nfeat*j];
				Z[i+nfeat*k] = Z[i+nfeat*k] + P[k+nclass*j];
			}
		}
	}
	// normalize C by sum P to have mean and std of each cluster
	for (i = 0; i < nfeat; i++) {
		Sbar = 0;
		for (k = 0; k < nclass; k++) {	
			C[i+nfeat*k] = C[i+nfeat*k] / Z[i+nfeat*k];
			S[i+nfeat*k] = maxOf(0.5, sqrt(S[i+nfeat*k] / Z[i+nfeat*k] - C[i+nfeat*k]*C[i+nfeat*k]));
			Sbar += S[i+nfeat*k]; 
		}
		Sbar = Sbar / nclass;
		for (k = 0; k < nclass; k++) {
			S[i+nfeat*k] = 0.75* S[i+nfeat*k] + 0.25 * Sbar;
		}
	}
}

function gmmUpdateImage(id,C,S,nclass) {
	// update the visualization of the classes
	if (id<0) {
		selectImage(id);
		nfeat=2;
		w = getWidth();
		h = getHeight();
		npts = w*h;
		X = newArray(nfeat*w*h);
		for (k=0;k<nclass;k++) {
			for (y = 0; y < h; y++) {
				for (x = 0; x < w; x++) {
					j = x+w*y;
					X[nfeat*j] = -4+15*x/w;
					X[1+nfeat*j] = (-4+15*(h-y)/h);
				}
			}
		}
		P = newArray(nclass*w*h);
		gmmUpdateProbabilities(X,C,S,P,npts,nfeat,nclass,1);
		for (k=0;k<nclass;k++) {
			Stack.setChannel(k+1);
			for (y = 0; y < h; y++) {			
				for (x = 0; x < w; x++) {
					j = x+w*y;
					setPixel(x,y,P[k+nclass*j]);
				}
			}
		}
	}
}

function gmmComputeStatsByClass(X,K,npts,nfeat,nclass) {
	/** Compute count and mean of the feature X for each class
	 * 
	 * X : features vector
	 * K : class index in [0,nclass[
	 * 
	 * return count, and mean of each features for each class (col major)
	 * class1 : count, feat1, feat2, feat3
	 * class2 : count, feat1, feat2, feat3
	 */	 
	cnt = newArray((nfeat+1)*nclass);
	for (j = 0; j < npts; j++) {
		k = K[j];
		cnt[(nfeat+1)*k] = cnt[(nfeat+1)*k] + 1;
		for (i = 0; i < nfeat; i++) {
			cnt[1+i+(nfeat+1)*k] = cnt[1+i+(nfeat+1)*k] + X[i+nfeat*j];
		}
	}
	// normalize the sum by the count store in (0,k)
	for (k = 0; k < nclass; k++) {
		for (i = 0; i < nfeat; i++) {
			cnt[1+i+(nfeat+1)*k] = cnt[1+i+(nfeat+1)*k] / cnt[(nfeat+1)*k];
		}
	}
	return cnt;
}

function gmmSaveStatsToTable(tblname,filename,group,stats,features,classes) {
	/* Save stats into a table
	 * Parameters
	 * tblname : name of the table
	 * group : goup (or condition)
	 * stats : array of stats
	 * features :array of features
	 * classes : array with classes names
	 */
	if (!isOpen(tblname)) {
		Table.create(tblname);
	}
	selectWindow(tblname); 
	nfeat = features.length;
	row = Table.size;
	for (k = 0; k < classes.length; k++) {
		Table.set("File", row, filename);
		Table.set("Group", row, group);
		Table.set("Class", row, classes[k]);
		Table.set("Count", row, stats[(nfeat+1)*k]);
		for (i = 0; i < nfeat; i++) {
			Table.set(features[i], row, stats[i+1+(nfeat+1)*k]);
		}
		row++;
	}
	Table.update;
}

function gmmPrintStats(stats, features, classes)  {
	nfeat = features.length;
	for (k = 0; k < classes.length; k++){
		str = classes[k]+" : n:"+ stats[(nfeat+1)*k];
		for (i = 0; i < nfeat; i++) {
			str += " "+features[i]+":"+stats[i+1+(nfeat+1)*k];
		}
		print(str);
	}
}

function filterROI(X,npts,nfeat,min,max) {
	/* Filter out ROI based on features within min x max hypercube
	 * 
	 */
	I = newArray(npts); // 1 if set to delete
	n = 0; // number of ROI to keep
	for (j = 0; j < npts; j++) {
		for (i = 0; i < nfeat; i++) {
			if (X[i+nfeat*j] < min(i) || X[i+nfeat*j] > max(i)) {
				I[j] = 1;
			}
		}
		if (I[j]==0) {
			n++;
		}
	}
	Y = newArray(n*nfeat);
	jj = 0;
	for (j = npts; j > 0; j--) {
		if (I[j] == 1) {
			roiManager("select", j);
			roiManager("delete");
		} else {
			for (i = 0; i < nfeat; i++) {
				Y[i+nfeat*jj] = X[i+nfeat*j];
			} 
			jj++;
		}
	}
	return Y;
}