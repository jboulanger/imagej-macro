//setBatchMode("hide");
print("\\Clear");
Overlay.remove();

f = getSignalFit(0,500,100,2);
Plot.create("Title", "X-axis Label", "Y-axis Label",Array.getSequence(f.length),f);
/*
x = getValue("X");
y = getValue("Y");
print("x:"+x+", y:"+y);




position = newArray(x,y,100,2.5,1);
T =  nSlices;
I = newArray(T);
J = newArray(T);
B = newArray(T);
for (frame = 1; frame <= T; frame++) {
	setSlice(frame);
	position = measureSpotIntensity(position);
	makePoint(position[0], position[1]);
	Overlay.addSelection;
	Overlay.setPosition(frame);
	Overlay.show;
	I[frame-1] = position[4];
	makeRectangle(x-2, y-2, 5, 5);
	J[frame-1] = getValue("Mean")-100;
	B[frame-1] = position[2];
}
Array.print(B);
Array.print(I);

Plot.create("Title", "X-axis Label", "Y-axis Label");
Plot.setColor("green");
Plot.add("line",Array.getSequence(T), I);
Plot.setColor("red");
Plot.add("line",Array.getSequence(T), J);
Plot.setColor("blue");
Plot.add("line",Array.getSequence(T), B);
Plot.setLimitsToFit();
Plot.update();
*/

//setBatchMode("exit and display");

function getSignalFit(idx,n,bg,sigma) {
	roiManager("select", idx);
	x = getValue("X");
	y = getValue("Y");	
	position = newArray(x,y,bg,sigma,1);
	f = newArray(n);	
	for (i = 0; i < n; i++) {
		Stack.setFrame(i + 1);
		position = estimateSpotParameters(position);
		f[i] = position[4];
		print("Position " + position[0] + "," + position[1]);
		makePoint(position[0], position[1]);
		Overlay.addSelection;
		Overlay.setPosition(i+1);
		Overlay.show;		
	}	
	return f;
}

function estimateSpotParameters(position) {
	x = position[0];
	y = position[1];
	bg = position[2];
	sigma = position[3];
	A = position[4];
	h = Math.ceil(3 * sigma);
	N = 10;
	for (i = 1; i < N; i++) {
		swx = 0;
		swy = 0;
		sw = 0;
		swb = 0;
		sb = 0;
		n = 0;
		for (yi = round(y - h); yi <= round(y + h); yi++) {
			for (xi = round(x - h); xi <= round(x + h); xi++) {
				dx = (xi - x) / sigma;
				dy = (yi - y) / sigma;
				w0 = exp(-0.5*(dx*dx+dy*dy));
				I = getPixel(xi, yi);
				w = w0 * maxOf(0, I - bg);
				swx += w * xi;
				swy += w * yi;
				sw += w;
				n += w0;
				sb += (1 - w0) * ( I - A * w0);
				swb += (1 - w0);
			}
		}
		if (sw / n > 10)  {
			x = swx / sw;
			y = swy / sw;
		}
		if (swb > 1) {
			bg = sb / swb;
		}
		A = sw/n;
	}
	return newArray(x,y,bg,sigma,A);
}
