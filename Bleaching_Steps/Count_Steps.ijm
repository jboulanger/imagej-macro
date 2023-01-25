//@File(label="Input file", description="csv or tiff+zip file") input
//@Float(label="Sensor gain", value=0.25) gain0
//@Float(label="Sensor offset", value=100) offset0
//@Float(label="Sensor readout", value=1) sigma0
//@Integer(value=10) iterations

/*
 * Counting bleaching steps
 *
 * Input:
 * Input can be a tif + zip file with ROIs or a csv file with signals.
 * 
 * If an image is opened and the input is not ending in csv or tif, process the current image and rois.
 * If there is only one ROI a grpahis produced at the end.
 * 
 * Output:
 * It outputs two tables
 * 
 * [ smoothing.csv ]
 *  For each roi the table has a pair of colum signal and the corresponding approximation.
 * 
 * [ step-statistics.csv ]
 * This table summarize the step statistics
 * - Positive step counts : number of positive steps > threshold
 * - Positive step total: number of positive steps > threshold
 * - Positive step mean: mean of positive steps > threshold
 * - Positive step median: median of positive steps > threshold
 * - Negative step counts : number of negative steps > threshold
 * - Negative step total: number of negative steps > threshold
 * - Negative step mean: mean of negative steps > threshold
 * - Negative step median: median of negative steps > threshold
 * - Negative step corrected count: ratio Negative step total / Negative step median
 * - Negative step rate : rate of the negative steps events
 * - Noise std : noise estimated from the signal after variance stabilization (should be 1)
 
 * Jerome Boulanger for Alex Fellow 2023
 */

start_time = getTime();
setBatchMode("hide");


if (matches(input,".*test")) {
	print("*** Testing mode ***");
	test_noisestd();
	test_gat();
	test_process();
	exit;
}

if (endsWith(input, ".csv")) {
	process_csv(input);
} else if (endsWith(input, ".tif")) {
	process_tif_and_roi(input);
} else {
	process_active_image();
}

setBatchMode("show");
print("Elapsed time " + (getTime() - start_time)/1000 + " seconds.");

function process_active_image() {

	print("Processing current image " + getTitle());
	Stack.getDimensions(width, height, channels, slices, frames);

	tbl2 = "smoothing.csv";
	Table.create(tbl2);

	tbl3 = "step-statistics.csv";
	Table.create(tbl3);

	N = roiManager("count");

	showStatus("Processing ROIs");
	for (i = 0; i < N; i++) {
		showProgress(i, N-1);
		roiManager("select", i);
		roiname = Roi.getName;
		f0 = getSignal(i, frames);
		f1 = gat(f0, gain0, offset0, sigma0);
		sigma = estimateNoiseStd(f1);
		u1 = smoothSteps(f1, sigma, iterations);
		u0 = igat(u1,  gain0, offset0, sigma0);		
		saveResult(tbl2, roiname, f0, u0);		
		summarizeSteps(tbl3, i, roiname, u0, sigma);		
	}
		
	if (N == 1) {
		print("Make a figure");
		makeFigure(f0, u0);
	}
}

function process_csv(input) {
	/* Process a csv file with extracted 1d signals
	 *  input (string) : path to the csv file
	 */
	print("Processing csv file " + input);
	Table.open(input);
	tbl1 = Table.title;
	headings = split(Table.headings,"\t");

	N = headings.length;
	N = 10;

	tbl2 = replace(File.getName(input), ".csv", "-smoothing.csv");
	Table.create(tbl2);
	
	tbl3 = replace(File.getName(input), ".csv", "-step-statistics.csv");
	Table.create(tbl3);

	showStatus("Processing ROIs");
	for (i = 0; i < N; i++) {
		showProgress(i, N-1);
		selectWindow(tbl1);
		print(headings[i]);
		f0 = Table.getColumn(headings[i]);
		f1 = gat(f0, gain0, offset0, sigma0);
		sigma = estimateNoiseStd(f1);
		u1 = smoothSteps(f1, sigma, iterations);
		u0 = igat(u1,  gain0, offset0, sigma0);		
		saveResult(tbl2, headings[i], f1, u1);
		summarizeSteps(tbl3, i, headings[i], u0, sigma);
	}
		
}

function process_tif_and_roi(input) {
	/* Process a tiff file and associated list of ROI
	 *  input (string): path to the file
	 */
	print("Processing tif file " + input);
	open(input);
	Stack.getDimensions(width, height, channels, slices, frames);

	roifile =  replace(input,".tif","_ROI.zip")
	if (!File.exists(roifile)) {
		error("Cound not find ROIs file " + roifile);
	}
	roiManager("Open", roifile);
	N = roiManager("count");

	tbl2 = replace(File.getName(input), ".csv", "-smoothing.csv");
	Table.create(tbl2);

	tbl3 = replace(File.getName(input), ".csv", "-step-statistics.csv");
	Table.create(tbl3);

	showStatus("Processing ROIs");
	for (i = 0; i < N; i++) {
		showProgress(i, N-1);
		roiManager("select", i);
		roiname = Roi.getName;
		f0 = getSignal(i, frames);
		f1 = gat(f0, gain0, offset0, sigma0);
		sigma = estimateNoiseStd(f1);
		u1 = smoothSteps(f1, sigma, iterations);
		u0 = igat(u1,  gain0, offset0, sigma0);		
		saveResult(tbl2, roiname, f0, u0);
		summarizeSteps(tbl3,i, roiname, u0, sigma);
	}
}

function saveResult(tbl, roiname, f, u) {
	/* Save the filtered signal and the original one in a table
	 * tbl (string) : table name
	 * f (array) : original signal
	 * u (array) :filtered signal
	 */
	selectWindow(tbl);
	n = f.length;
	for (i = 0; i < n; i++) {
		Table.set(roiname+"-signal", i, f[i]);
		Table.set(roiname+"-approx", i, u[i]);
	}	
	Table.update;
}

function summarizeSteps(tbl, idx, roiname, f, sigma) {
	/* Sumarize step counting statistics in a table
	 *  tbl (string): table name
	 *  idx (number): index of the roi
	 *  roiname (string) : name of the roi
	 *  f (array) : 1d filtered signal
	 */
	selectWindow(tbl);
	threshold = 0.125;
	n = f.length;
	positive_steps_amplitude = newArray(n);
	positive_steps_indices = newArray(n);
	positive_steps_n = 0;
	positive_steps_total = 0;
	
	negative_steps_amplitude = newArray(n);
	negative_steps_indices = newArray(n);
	negative_steps_n = 0;
	negative_steps_total = 0;
	
	for (i = 0; i < n-1; i++) {
		d = f[i+1] - f[i];
		if (d > threshold) {
			positive_steps_amplitude[positive_steps_n] = d;
			positive_steps_indices[positive_steps_n] = i;
			positive_steps_n++;
			positive_steps_total += d;
		} else if (d < - threshold){
			negative_steps_amplitude[negative_steps_n] = -d;
			negative_steps_indices[negative_steps_n] = i;
			negative_steps_n++;
			negative_steps_total -= d;
		}
	}
	
	positive_steps_amplitude = Array.trim(positive_steps_amplitude, positive_steps_n);
	positive_steps_indices = Array.trim(positive_steps_indices, positive_steps_n);
	positive_step_size = median(positive_steps_amplitude);
	
	negative_steps_amplitude = Array.trim(negative_steps_amplitude, negative_steps_n);
	negative_steps_indices = Array.trim(negative_steps_indices, negative_steps_n);
	negative_step_size = median(negative_steps_amplitude);
	negative_step_count = round(negative_steps_total / negative_step_size);
	negative_step_rate = compute_rate(negative_steps_indices);
	
	selectWindow(tbl);
	Table.set("ROI", idx, roiname);	
	Table.set("Positive steps counts", idx, positive_steps_n);
	Table.set("Positive steps total", idx, positive_steps_total);
	Table.set("Positive steps mean", idx, positive_steps_total / (positive_steps_n));
	Table.set("Positive steps median", idx, positive_step_size);
	
	Table.set("Negative steps counts", idx, negative_steps_n);
	Table.set("Negative steps total", idx, negative_steps_total);
	Table.set("Negative steps mean", idx, negative_steps_total / (negative_steps_n));
	Table.set("Negative steps median", idx, negative_step_size);
	Table.set("Negative steps corrected count", idx, negative_step_count);
	Table.set("Negative steps rate", idx, negative_step_rate);
	
	Table.set("Noise std dev", idx, sigma);
	
	Table.update;
}

function compute_rate(indices) {
	/* Compute the bleaching rate from the bleaching step indices
	 *  indices (array): list of the steps
	 *  Returns (number) the rate at which the steps
	 */
	n = indices.length;
	if (n==0) {
		return 0;
	}
	R = newArray(n-1);
	for (i = 0; i < n - 1; i++) {
		R[i] = (n - i) * (indices[i+1] - indices[i]);
	}
	rate = 1 / median(R);
	return rate;
}


function makeFigure(f,u) {
	/* Make a graph with the signal and its approximation
	 *  f (array) : signal
	 *  u (array) : approximation
	 */
	x = Array.getSequence(f.length);
	Plot.create("Bleaching steps", "Frame", "Intensity");
	Plot.setColor("black");
	Plot.setLineWidth(1);
	Plot.add("line", x, f);
	
	Plot.setColor("red");
	Plot.setLineWidth(2);
	Plot.add("line", x, u);
	Plot.update();
}

function getSignal(idx, n) {	
	/* Extract the intensity over time for the cirrent ROI
	 *  idx (number) : index of the ROI
	 *  n (number) : number of frames
	 */
	f = newArray(n);
	for (i = 0; i < n; i++) {
		roiManager("select", idx);
		Stack.setFrame(i + 1);
		f[i] = getValue("Mean");
	}	
	return f;
}

function smoothSteps(f,sigma,niter) {
	/* Smooth a 1d signal in a step function
	 *  f (array) : original signal
	 *  sigma (number) : noise standard deviation
	 *  niter (number) : number of iterations
	 * Returns (array) filtered signal
	 */
	n = f.length;
	u0 = Array.copy(f);
	u1 = newArray(n);
	v0 = newArray(n);
	v1 = newArray(n);

	for (i = 0; i < n; i++) {
		v0[i] = sigma * sigma;
	}

	for (iter = 0; iter < niter; iter++) {

		h = pow(2.0, iter);

		for (i = 0; i < n; i++) {

			j0 = maxOf(0, i - h);
			j1 = minOf(n - 1, i + h);

			swu = 0;
			sw = 0;

			for (j = j0; j <= j1; j++) {

				d = (u0[j] - u0[i]) /(6 * sqrt(v0[i]));
				if (abs(d) < 1) {
					sw += 1;
					swu += f[j];
				}
			}

			if (sw > 0) {
				u1[i] = swu / sw;
				v1[i] = sigma * sw / (sw*sw);
			}
		}

		u0 = Array.copy(u1);
		v0 = Array.copy(v1);
	}
	return u1;
}

function median(x) {
	/* Median
	 *  x (array): list of values
	 * Return the median of the array x
	 */
	y = Array.copy(x);
	Array.sort(y);
	n = y.length;
	if (n==0) {
		return 0;
	}
	if (n%2==0) {
		return 0.5*(y[n/2-1] + y[n/2]);
	} else {
		return y[(n-1)/2];
	}
}

function estimateNoiseStd(f) {
	/* Estimate noise standard deviation
	 *  f (array) : 1d signal
	 * Returns the estimated noise standard deviation
	 */
	n = f.length-1;
	delta = newArray(n);	
	for (i = 0; i < n; i++) {
		delta[i] = f[i+1] - f[i];	
	}
	m = median(delta);
	for (i = 0; i < n; i++) {
		delta[i] = abs(delta[i] - m);
	}
	m = 1.48 * median(delta) / sqrt(2);
	return m;
}

function test_noisestd() {
	/*Test the noise std estimation*/
	n = 1000;
	f = newArray(n);
	for (i = 0; i < n; i++) {
		f[i] = random("gaussian");
	}	
	sigma = estimateNoiseStd(f);
	if (abs(sigma-1) < 0.1) {
		print("Noise std estimation test ok");
	} else {
		print("Noise std estimation test failed sigma:"+sigma);
	}	
}

function test_gat() {
	/* Test the gat functions */
	n = 1000;
	f = newArray(n);
	g = 0.2;
	m = 100;
	s = 1;
	for (i = 0; i < n; i++) {
		v = 10 + i;
		f[i] = g * (v + sqrt(v) * random("gaussian")) + m + s*random("gaussian");
	}	
	u = igat(gat(f,g,m,s),g,m,s);
	ok = true;
	for (i = 0; i < f.length; i++) {
		ok = ok && (abs(f[i] - u[i]) < 0.01);
	}
	if (ok) print("GAT test ok");
	else print("GAT test failed");
	return ok;
}

function test_process() {
	N0 = 10;
	T = 3000;
	rate = 1e-3;
	gain0 = 0.2;
	offset0 = 100;
	sigma0 = 1;
	alpha = 100;
	
	// simulation of the step function
	u = newArray(T);	
	N = N0;
	for (i = 0; i < T; i++) {
		u[i] = N;		
		for (n = 0; n < N; n++) {			
			if (random < rate) {
				N--;
			}
		}
	}
	
	// simulation of the noise
	f0 = newArray(T);
	ubar = newArray(T);
	for (i = 0; i < T; i++) {
		v = alpha * u[i];
		f0[i] = gain0 * (v + sqrt(v) * random("gaussian")) + offset0 + sigma0 * random("gaussian");
		ubar[i] = gain0 * v + offset0;
	}	
	
	f1 = gat(f0, gain0, offset0, sigma0);
	sigma = estimateNoiseStd(f1);
	u1 = smoothSteps(f1, sigma, iterations);
	u0 = igat(u1,  gain0, offset0, sigma0);
	
	Table.create("test");
	summarizeSteps("test", 0, "test", u0, sigma);
	
	Plot.create("Bleaching steps test", "Time [frames]", "Intensity");
	x = Array.getSequence(T);
	Plot.setColor("#AAAAAA");
	Plot.add("line", x, f0);
	Plot.setColor("red");
	Plot.add("line", x, ubar);
	Plot.setColor("green");
	Plot.add("line", x, u0);
	
	selectWindow("test");
	Nestimated = Table.get("Negative steps counts", 0);	
	if (Nestimated == N0 - N) {
		print("Process test ok N:"+ Nestimated + "/" + N0-N);
	} else {
		print("Process test failed N:"+ Nestimated + "/" + N0-N);
	}
	
}

function gat(f,g,m,s) {
	/* Generalized Anscombe transform
	 *  f (array) : 1d signal
	 *  g (number) : gain
	 *  m (number) : offset
	 *  s (number) : standard deviation
	 * Return an array
	 */
	u = Array.copy(f);
	n = f.length;
	for (i = 0; i < n; i++) {
		a = 3 / 8 * g * g + s * s - g * m;
		u[i] = 2 / g * sqrt(maxOf(0, g * f[i] + a));
	}
	return u;
}

function igat(f,g,m,s) {
	/* Inverse generalized Anscombe transform
	 *  f (array) : 1d signal
	 *  g (number) : gain
	 *  m (number) : offset
	 *  s (number) : standard deviation
	 * Return an array
	 */
	u = Array.copy(f);
	n = f.length;
	for (i = 0; i < n; i++) {
		a = 3 / 8 * g * g + s * s  - g * m;
		u[i] = (g * g / 4 * f[i] *  f[i] - a ) / g;
	}

	return u;
}
