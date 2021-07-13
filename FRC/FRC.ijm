//@ImagePlus(label="first image") im1
//@ImagePlus(label="second image") im2
//@Boolean(label="add noise", value=false) addnoise
//@Boolean(label="use estimated noise", value=false) estnoise
//@Float(label="specified noise level", value=0.0) noise_level

/*
	Compute a Fourier Ring Correlatin curve
	Estimate a resolution using the 1/7 rule

	Jerome Boulanger (2008)
*/

macro "FRC" {	
	if ((im1 == im2) && !addnoise) {
		print("The FRC computed on a single image will not be informative unless you add noise");
	}
	setBatchMode(true);		
	selectImage(im1);
	name1 = getTitle;		
	run("Duplicate...", "title=u");
	run("32-bit");
	sigma1 = estimateNoiseLevel();	
	getPixelSize(unit, pw, ph, pd);
	List.setMeasurements;
	run("Subtract...", "value="+List.getValue("Mean")+" slice");	
	if (addnoise) {
		if (estnoise) {
			sigma = sigma1;
		} else { 
			sigma = noise_level;
		}
		run("Add Specified Noise...", "standard="+sigma);
	}
	// apply a window
	run("Macro...", "code=v=v*exp(-d*d/(0.1*w*w))");
	run("FFT Options...", "complex do");
	run("Stack to Images");	
	selectWindow("Real"); rename("a");
	selectWindow("Imaginary"); rename("b");	
	U = getImageID();
	
	selectImage(im2);	
	name2 = getTitle;
	run("Duplicate...", "title=v");
	run("32-bit");	
	sigma2 = estimateNoiseLevel();
	List.setMeasurements;	
	run("Subtract...", "value="+List.getValue("Mean")+" slice");
	if (addnoise) {
		if (estnoise) {
			sigma = sigma2;			
		} else { 
			sigma = noise_level;
		}
		run("Add Specified Noise...", "standard="+sigma);
	}	
	// apply a window
	run("Macro...", "code=v=v*exp(-d*d/(0.05*w*w))");
	run("FFT Options...", "complex do");
	run("Stack to Images");
	selectWindow("Real"); rename("c");
	selectWindow("Imaginary"); rename("d");
	V = getImageID();
	
	
	// compute sum(ac+bd ) / sqrt(sum(aa+bb)*sum(cc+dd))
	imageCalculator("Multiply create", "a","c"); rename("ac");
	imageCalculator("Multiply create", "b","d"); rename("bd");
	imageCalculator("Add", "ac", "bd");rename("ac+bd");
	selectWindow("a"); run("Square");
	selectWindow("b"); run("Square");
	selectWindow("c"); run("Square");
	selectWindow("d"); run("Square");
	imageCalculator("Add", "a", "b");rename("aa+bb");
	imageCalculator("Add", "c", "d");rename("cc+dd");
	selectWindow("ac+bd");
	cx = getWidth() / 2;
	cy = getHeight() / 2;	
	fmax = getWidth() / 2;
	N = getWidth() / 2;
	M = "Mean";
	radius = generateSequence(N, 0, fmax);
	UV = radialMeasure(cx,cy,radius,M);
	selectWindow("aa+bb");
	U2 = radialMeasure(cx,cy,radius,M);
	selectWindow("cc+dd");
	V2 = radialMeasure(cx,cy,radius,M);	
	UV[0] = 1;	
	Array.getStatistics(U2, mini, maxi, mean, stdDev);
	for (k = 1; k < UV.length; k++) {
		d = U2[k] * V2[k];
		if (d > 0) {
			UV[k] = (UV[k]) / sqrt(d);		
		} else {
			UV[k] = 0;
		}
		U2[k] = log(0.0001+U2[k]);
	}
	U2[0] = log(0.0001+U2[0]);
	Array.getStatistics(U2, mini, maxi, mean, stdDev);
	for (k = 0; k < UV.length; k++) {
		U2[k] = (U2[k]-mini) / (maxi-mini);
	}	
	UVS = smoothCurve(UV,5);
	Plot.create("FRC of "+name1+"/"+name2, "radius", "Amplitude");
	Plot.setColor("#34495e");
	Plot.add("line",radius, UV);
	Plot.setColor("#e74c3c");
	Plot.add("line",radius, UVS);
	Plot.setColor("#2ecc71");
	Plot.add("line",radius, U2);
	Plot.setColor("#3498db");
	Plot.drawLine(0, 1/7, N, 1/7);
	Plot.setLegend("FRC\nSmooth FRC\nPower Spectrum\n1/7\n");
	closeWindow("ac+bd");
	closeWindow("aa+bb");
	closeWindow("cc+dd");
	closeWindow("b");
	closeWindow("d");
	closeWindow("bd");	
	closeWindow("u");
	closeWindow("v");
	
	fcut = intersecCurve(UVS, 1/7) / (2 * N * pw);	
	i=nResults;
	setResult("Image 1",i,name1);
	setResult("Image 2",i,name2);
	setResult("Noise 1",i,sigma1);
	setResult("Noise 2",i,sigma2);	
	setResult("resolution (1/7)",i, "" +1 / fcut+" " + unit);
	
	setBatchMode(false);
}


macro "Radial Average" {
	if (nImages()==0) {
		newImage("Untitled", "32-bit black", 256, 256, 1);		
		r=64;
		makeOval(128-r,128-r,2*r,2*r);
		setColor(128,128,128);
		fill();
	}
	x = generateSequence(getWidth()/8,1,getWidth()/2);	
	y = radialMeasure(getWidth()/2,getHeight()/2, x, "Mean");	
	Plot.create("a","x","y",x,y);
}

/*
 * Compute the radial average of the current image 
 */
function radialMeasure(_cx,_cy, _radius,_measure) {	
	getPixelSize(_unit, _pw, _ph, _pd);	
	_avg = newArray(_radius.length);
	makeOval(_cx-_radius[1], _cy-_radius[1], 2*_radius[1], 2*_radius[1]);
	roiManager("Add");
	List.setMeasurements;
	_avg[0] = List.getValue(_measure);
	for (_k = 1; _k < _radius.length; _k++) {
		_dr = abs(_radius[_k] - _radius[_k-1]);		
		makeOval(_cx-_radius[_k], _cy-_radius[_k], 2*_radius[_k], 2*_radius[_k]);		
		run("Make Band...", "band="+(_dr*_pw));		
		List.setMeasurements;
		_avg[_k] = List.getValue(_measure);
	}
	run("Select None");	
	return _avg;
}

function generateSequence(_N,_start,_end) {
	_s = newArray(_N);
	for (_n = 0; _n < _N; _n++) {
		_s[_n] = _start + _n * (_end-_start) / (_N-1);
	}
	return _s;
}

function closeWindow(_str) {
	selectWindow(_str);close();
};

function smoothCurve(X,h) {
	Y = newArray(X.length);	
	for (i=0;i<X.length;i++) {
		acc = 0;
		n = 0;
		for (j=maxOf(0,i-3*h);j<minOf(X.length,i+3*h+1);j++) {	
			w = exp(-0.5*(i-j)*(i-j)/(h*h));
			acc = acc + w * X[j]; 
			n = n + w;
		}
		Y[i] = acc / n;
	}
	return Y;
}

function intersecCurve(X,value) {
	k = 0;
	while (k < X.length-1 && X[k] > value) {
		k++;
	}
	return k;
}

function estimateNoiseLevel() {
	run("Duplicate...","title=noise");
	run("Convolve...", "text1=[0.0 0.22361 0.0\n0.22361 -0.89443 0.22361\n0.0 0.22361 0.0\n]");	
	getStatistics(area, mean, min, max, std);
	close();
	return std;
}
