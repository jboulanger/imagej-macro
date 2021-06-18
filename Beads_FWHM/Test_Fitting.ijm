print("\\Clear");
newImage("Untitled", "32-bit black", 256, 256, 7);
setSlice(4);
setPixel(128,128,255);
run("Gaussian Blur 3D...", "x=2 y=2 z=2");
profile = circularAverage(128,128,12,9);
radius = getRadiusValues(profile.length);
Fit.doFit("Gaussian (no offset)", radius, profile, newArray(1,0,2));	
		a = Fit.p(0);
		b = Fit.p(1);
		c = Fit.p(2);
		print(c);
vals = circularFit(128,128,12,9);
Array.getStatistics(vals, min, max, mean, stdDev);
print("mean:"+mean);

function circularFit(x0,y0,N,L) {
	getPixelSize(unit, pw, ph, pd);
	makeOval(x0-L+0.5,y0-L+0.5,2*L,2*L);	
	run("Add Selection...");
	vals = newArray(N);
	for (i = 0; i < N; i++) {
		makeLine(x0-L*cos(PI*i/N), y0-L*sin(PI*i/N),x0+L*cos(PI*i/N), y0+L*sin(PI*i/N));		
		run("Add Selection...");
		profile = getProfile();
		radius = getRadiusValues(profile.length);
		Fit.doFit("Gaussian (no offset)", radius, profile, newArray(1,0,2));		
		vals[i]= Fit.p(2);
		print(vals[i]);
	}
	run("Select None");
	return vals;
}

function circularAverage(x0,y0,N,L) {
	getPixelSize(unit, pw, ph, pd);
	makeOval(x0-L+0.5,y0-L+0.5,2*L,2*L);	
	run("Add Selection...");
	for (i = 0; i < N; i++) {
		makeLine(x0-L*cos(PI*i/N), y0-L*sin(PI*i/N),x0+L*cos(PI*i/N), y0+L*sin(PI*i/N));		
		run("Add Selection...");
		profile = getProfile();
		if (i==0) {
			avg = newArray(profile.length);
		}
		for (k = 0; k < profile.length && k < avg.length; k++) {
			avg[k] = avg[k] + profile[k];
		}
	}
	Array.getStatistics(avg, ymin, ymax);
	for (k = 0; k < avg.length; k++) {
		avg[k] = (avg[k] - ymin) / (ymax-ymin);
	}
	run("Select None");
	Overlay.show();
	return avg;
}

function getRadiusValues(N) {	
	radius = newArray(N);			
	for (k = 0; k < N; k++) {
		radius[k] = (k-radius.length/2);			
	}
	return radius;
}