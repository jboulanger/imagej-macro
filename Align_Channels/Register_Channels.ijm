#@String(label="mode",choices={"Estimate","Correct"}) mode
/* Register channels using Motion estimation
 * 
 */

if (nImages==0) {
	print("\\Clear");
	print("Demo");
	run("Boats");
	run("Copy");
	run("Add Slice");
	run("Paste");
	print("True shift: 5,10 ad -5,-7");
	run("Translate...", "x="+5+" y="+10+" interpolation=Bicubic slice");
	run("Add Slice");
	run("Paste");
	run("Translate...", "x="+-3+" y="+-7+" interpolation=Bicubic slice");
	run("Stack to Hyperstack...", "order=xyczt(default) channels=3 slices=1 frames=1 display=Composite");
	run("Add Noise");	
}	

start = getTime();
setBatchMode(true);
if (matches(mode, "Estimate")) {	
	align_channel();
} else {
	correct();
}
setBatchMode(false);
print("Elapsed time "+ (getTime()-start)/1000 + "s");


function correct() {
	// Correct the channel drift from the stack using the value in the table
	Stack.getDimensions(width, height, channels, slices, frames);
	for (channel = 2; channel <= channels; channel++){
		showProgress(channel,channels);
		vx = getResult("X", channel - 1);
		vy = getResult("Y", channel - 1);
			
		if (vx*vx+vy*vy>0.0001) {
			for (slice = 1; slice <= slices; slice++) {
				for (frame = 1; frame <= frames; frame++) {					
					Stack.setPosition(channel, slice, frame);
					run("Translate...", "x="+-vx+" y="+-vy+" interpolation=Bicubic slice");
				}
			}
		}
	}
}

function align_channel() {
	id = getImageID();
	Stack.getDimensions(width, height, channels, slices, frames);
	for (k=0;k<channels;k++){setResult("Channel", k, k+1);setResult("X",k,0);setResult("Y",k,0);}
	for (i =  1; i <= 1; i++) {
		for (j = i+1; j <= channels; j++) if(i!=j) {
			selectImage(id);
			run("Duplicate...", "title=A duplicate channels="+i);
			id1 = getImageID();
			selectImage(id);
			run("Duplicate...", "title=B duplicate channels="+j);
			id2 = getImageID();
			shift = estimateShift(id1,id2);
			print("shift "+i+"-"+j+": "+shift[0]+", "+shift[1]);
			setResult("Channel", j-1, j);
			setResult("X", j-1, getResult("X", j-1)+shift[0]);
			setResult("Y", j-1, getResult("Y", j-1)+shift[1]);
			selectImage(id);
			Stack.setPosition(j, 1, 1);
			run("Translate...", "x="+-shift[0]+" y="+-shift[1]+" interpolation=Bicubic slice");
			selectImage(id1);close();
			selectImage(id2);close();
		}
	}
}


function estimateShift(id1,id2) {
	sx = 0;
	sy = 0;
	n = 1;
	for (j = 0; j < 3; j++) {
		scale = pow(2,j-2);
		for (iter = 0; iter < 100 && n > 0.01; iter++) {
			selectImage(id2); run("Duplicate...", " "); id3 = getImageID();
			run("Translate...", "x="+-sx+" y="+-sy+" interpolation=Bilinear");
			shift = estimateShiftCore(id1,id3,scale);
			selectImage(id3); close();
			sx -= shift[0];
			sy -= shift[1];
			n = sqrt(shift[0]*shift[0]+shift[1]*shift[1]);
		}
	}
	return newArray(sx,sy);
}

function estimateShiftCore(id1,id2,scale) {
	selectImage(id1);
	w = round(getWidth()*scale);
	h = round(getHeight()*scale);
	run("Duplicate...", " ");
	run("32-bit");
	if (scale != 1) {
		run("Size...", "width="+w+" height="+h+" depth=1 constrain average interpolation=Bilinear");	
	}
	ix = getImageID();rename("ix");
	selectImage(ix);
	run("Duplicate...", " ");
	iy = getImageID();rename("iy");

	imageCalculator("Subtract create 32-bit", id2, id1);
	it = getImageID();rename("it");
	if (scale > 0) {
		run("Size...", "width="+w+" height="+h+" depth=1 constrain average interpolation=Bilinear");	
	}

	selectImage(ix);
	run("Convolve...", "text1=[0 0 0\n-1 0 1\n0 0 0\n]");
	
	selectImage(iy);
	run("Convolve...", "text1=[0 -1 0\n0 0 0\n0 1 0\n]");
	// compute the second order images ixx iyy 
	imageCalculator("Multiply create 32-bit", ix, it); ixt = getImageID();rename("ixt");
	imageCalculator("Multiply create 32-bit", iy, it); iyt = getImageID();rename("iyt");
	imageCalculator("Multiply create 32-bit", ix, iy); ixy = getImageID();rename("ixy");
	selectImage(ix); run("Square");
	selectImage(iy); run("Square");
	// create a mask
	imageCalculator("Multiply create 32-bit", id1, id2); mask = getImageID();rename("mask");
	run("Macro...", "code=v=(v>1) slice");
	ids = newArray(ix,iy,it,ixy,ixt,iyt);
	for (i = 0; i < ids.length; i++) {imageCalculator("Multiply", ids[i], mask);}
	// measure the mean intensity in each image
	getDimensions(w, h, channels, slices, frames);
	band = 2;
	selectImage(ix);  makeRectangle(band, band, w-2*band, h-2*band); a = getValue("Mean");
	selectImage(ixy); makeRectangle(band, band, w-2*band, h-2*band); b = getValue("Mean");
	selectImage(iy);  makeRectangle(band, band, w-2*band, h-2*band); c = getValue("Mean");
	selectImage(ixt); makeRectangle(band, band, w-2*band, h-2*band); d = getValue("Mean");
	selectImage(iyt); makeRectangle(band, band, w-2*band, h-2*band); e = getValue("Mean");
	// solve [a,b;b,c] [vx;vy] = [d;e] 
	vx = 2*(c*d-b*e) / (a*c-b*b);
	vy = 2*(a*e-b*d) / (a*c-b*b);
	// clear the list of images
	for (i = 0; i < ids.length; i++) {selectImage(ids[i]);close();}
	selectImage(mask);close();
	return newArray(vx/scale,vy/scale);
}

