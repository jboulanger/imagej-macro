#@File (label="Input file") ipath_bigstack
#@String (label="Look up table", value="Red,Green,Blue") luts_str
#@Integer (label="Font Size", value=20) font_size

/*
 * MIP RGB Stack
 * 
 * Convert files by computing a MIP converting them to RGB and stacking them in one image
 * 
 * This is useful to get a quick overview of the data when combine with a montage
 * and using it in batch 
 * 
 * Jerome Boulanger
 */

var Width; // width of the stack
var Height; // height of the stack
var file_count; // number of files
var sid; // stack id

luts  = parseCSVString(luts_str);

if(file_count==0) {run("Close All");}

id = getFlattenedImage(ipath_bigstack, luts, font_size);

sid = addToStack(sid, id);

file_count++;

function addToStack(sid, id) {
	// add the 2D flattened image to a stack
	if (sid < 0) {			
		run("Copy");
		selectImage(sid);	
		setSlice(nSlices);
		run("Add Slice");
		run("Paste");		
		selectImage(id);
		close();	
	} else {
		rename("BigStack");	
		sid = getImageID();	
		getDimensions(Width, Height, dummy, dummy, dummy);
	}
	selectImage(sid);
	run("Select None");
	return sid;
}

function getFlattenedImage(path, luts,font_size) {
	// mip and color normalization, resize, overlay, flatten
	open(path);
	id0 = getImageID();	
	// mip
	run("Z Project...", "projection=[Max Intensity]");
	// color normalization
	Stack.getDimensions(dummy, dummy, channels, dummy, dummy);	
	NC = minOf(channels, luts.length);	
	for (c = 0 ; c < NC; c++) {	
		Stack.setChannel(c+1); 
		run(luts[c]); 
		autoContrast();
	}
	id1 = getImageID();
	// convert to RGB
	run("RGB Color");
	id2 = getImageID();
	// resize to the size of the first image in the batch
	if (Width != 0) {
		run("Size...", "width="+Width+" height="+Height+" depth=1 constrain average interpolation=Bilinear");	
	}	
	// add overlay
	setColor("white");
	setFont("Arial", 24);
	Overlay.drawString(File.getName(path), 5, 30);
	Overlay.show();
	// flatten
	run("Flatten");
	id3  = getImageID();
	selectImage(id0);close();
	selectImage(id1);close();
	selectImage(id2);close();
	return id3;
}

function autoContrast() {
	// adjust the contrast of the image based on basic statistics	
	getDimensions(w, h, dummy, dummy, dummy);
	npx = w*h;
	nbins = 255;
	getHistogram(values, counts, nbins);
	Array.getStatistics(counts, min, max, mean, stdDev);
	
	// find the mode of the histogram
	pos = Array.findMaxima(counts, max/10);	
	if (pos.length>0) {
		print("Using mode");
		vmin = values[pos[0]]; 
		sum = 0;
		for (i = 0; i < nbins && sum < 0.999 * npx; i++) {
			sum += counts[i];			
		}
		vmax = values[i-1];
	} else {
		getStatistics(area, mean, min, max, std, histogram);
		vmin = mean;
		vmax = mean+5*std;
	}	
	print("v=["+vmin+"-"+vmax+"]");
	setMinAndMax(vmin,vmax);
}

function parseCSVString(csv) {
	
	str = split(csv, ",");
	for (i = 0; i < str.length; i++) { 
		str[i] =  String.trim(str[i]);
	}
	return str;
}