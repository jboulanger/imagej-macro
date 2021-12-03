/*
 * Spot colocalization based on the distance between centroids
 * 
 * Jerome for Anna
 */
run("Select None");
rois = detectSpotsDOG(1, -1);


function objectColocalization(roi) {
	Stack.getDimensions(width, height, channels, slices, frames);
	n = roiManager("count");
	
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			ci = roi[3*i+2];
			cj = roi[3*j+2];
			if (c1!=cj) {
				dx = roi[3*i]-roi[3*j];
				dy = roi[3*i+1]-roi[3*j+1];
				d = dx*dx + dy*dy;
				if (d < 1) {
					
				}
			}
		}		
	}
}

function detectSpotsDOG(scale, logpfa) {
	// detect spots in all channels and returns an array x,y,c for each roi
	id0 = getImageID();
	run("Duplicate...", "duplicate");
	run("32-bit");
	id1 = getImageID();
	run("Gaussian Blur...", "sigma="+scale);
	run("Duplicate...", "duplicate");
	run("Gaussian Blur...", "sigma="+3*scale);
	id2 = getImageID();
	imageCalculator("Subtract 32-bit stack", id1,id2);
	selectImage(id1);
	Stack.getDimensions(width, height, channels, slices, frames);
	for (c = 1; c <= channels; c++) {
		Stack.setChannel(c);
		t = getValue("Mean")+getValue("StdDev");
		run("Macro...", "code=v=255*(v>"+t+") slice");
	}
	run("8-bit");
	run("Convert to Mask", "method=Default background=Light calculate");
	run("Watershed", "stack");
	roi = newArray();
	for (c = 1; c <= channels; c++) {
		Stack.setChannel(c);
		run("Select None");
		n0 = roiManager("count");
		run("Analyze Particles...", "size=1-Infinity pixel add slice");
		n1 = roiManager("count");
		tmp = newArray(3*(n1-n0));
		for (i = 0; i < n1-n0; i++) {
			tmp[3*i]   = getValue("X");
			tmp[3*i+1] = getValue("Y");
			tmp[3*i+2] = c;
		}
		roi = Array.concat(roi,tmp);
	}
	selectImage(id1);close();
	selectImage(id2);close();
	return roi;
}

function threshold(logpfa) {
	 // Define a threshold based on a probability of false alarm (<0.5)
	 st = laplaceparam();
	 pfa = pow(10,logpfa);
	 Array.print(st);
	 print("pfa " + pfa);
	 if (pfa > 0.5) {
	 	t = st[0] + st[1]*log(2*pfa);
	 } else {
	 	t = st[0] - st[1]*log(2*pfa);
	 }
	 run("Macro...", "code=v=255*(v<"+t+")");
}

function laplaceparam() {
	/** compute laplace distribution parameters 
	 *  Return center, scale and maximum value
	**/
	getHistogram(values, counts, 256);
	n = getWidth() * getHeight();
	c = 0;
	for (i = 0; i < counts.length && c <= n/2; i++) {
		c = c + counts[i];
	}
	m = values[i-1];
	b = 0;
	for (i = 0; i < counts.length; i++){
		b = b + abs(values[i] - m) * counts[i];
	}
	return newArray(m,b/n,values[values.length-1]);
}
