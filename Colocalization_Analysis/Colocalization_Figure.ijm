#@File(label="Input File") path
#@Integer(label="Channel 1",value=1) ch1
#@Float(label="Min 1",value=0) min1
#@Float(label="Max 1",value=255) max1
#@ String(choices={"Red","Magenta","Green","Cyan","Blue"}, style="dropdown") lut1
#@Integer(label="Channel 2",value=2) ch2
#@Float(label="Min 2",value=0) min2
#@Float(label="Max 2",value=255) max2
#@ String(choices={"Red","Magenta","Green","Cyan","Blue"}, style="dropdown") lut2
#@Boolean(label="Save PNG",value=True) savepng

/*
 * Colocalization figure
 *
 * Create a montage with two channels side by side
 * with each channels in gray scale and the composite in
 * a third column.
 *
 *
 * Jerome Boulanger for Kelly Nguyen (2024)
 */

function setActiveChannel(ch1, ch2) {
	/*
	 * Set pair of active channels in a composite image.
	 */
	str = "";
	Stack.getDimensions(width, height, channels, slices, frames);
	for (c = 1; c <= channels; c++) {
		if (c==ch1 || c==ch2) {
			str += "1";
		} else {
			str += "0";
		}
	}
	Stack.setActiveChannels(str);
}

function getPixelValues(img, channel) {
	/*
	 * Get pixel values in the 3D image
	 *
	 * Args:
	 *  img (int): id of the image
	 *  channel (int): index of the channel
	 *
	 * Returns
	 *  array of intensity values
	 *
	 */
	selectImage(img);
	Stack.getDimensions(width, height, channels, slices, frames);
	Stack.setDisplayMode("grayscale");
	values = newArray(width * height * slices);
	for (z = 1; z <= slices; z++) {
		Stack.setPosition(channel,z,1);
		for(y = 0; y < height; y++) {
			for(x = 0; x < width; x++) {
				values[x+width*(y+height*z)] = getPixel(x,y);
			}
		}
	}
	return values;
}


function PCC(ch1,ch2) {

	id = getImageID();
	
	selectImage(id);
	run("Duplicate...", "duplicate channels="+ch1);
	run("32-bit");
	run("Subtract...", "value="+getValue("Mean"));
	std1 = getValue("StdDev");
	id1 = getImageID();
	
	selectImage(id);
	run("Duplicate...", "duplicate channels="+ch2);
	run("32-bit");
	run("Subtract...", "value="+getValue("Mean"));
	std2 = getValue("StdDev");
	id2 = getImageID();
	
	imageCalculator("Multiply create 32-bit", id1, id2);
	pcc = getValue("Mean") / (std1*std2);
	print("PCC = " + pcc);
	close();
	selectImage(id);
	return pcc;
}


function hist2d(img, ch1, min1, max1, ch2, min2, max2) {
	/*
	 * Compute a 2D histogram
	 *
	 * Args:
	 * 	img (int): id of the image
	 * 	chs (array): indices of the channels
	 * 	
	 * Result:
	 *  int: index of the histogram
	 */
	 
	c1 = getPixelValues(img, ch1);
	c2 = getPixelValues(img, ch2);
	
	if (bitDepth()==8) {
		min1 = 0;
		min2 = 0;
		max1 = 255;
		max2 = 255;
	} else {
		Array.getStatistics(c1, min1, max1, mean, stdDev);
		Array.getStatistics(c2, min2, max2, mean, stdDev);
	}
	
	newImage("hist2d", "32-bit black", 255, 255, 1, 1, 1);
	for (i = 0; i < c1.length; i++) {
		u = round((c1[i] - min1) / (max1 - min1) * 255);
		v = round(255 - (c2[i] - min2) / (max2 - min2) * 255);
		val = getPixel(u, v) + 1;
		setPixel(u, v, val);
	}
	

	run("Square Root");
	run("Square Root");
	run("Enhance Contrast", "saturated=0.35");
	run("Gaussian Blur...", "sigma=1");
	run("Fire");
	return(getImageID());
}

function preprocess(img) {
	/*
	 * Subtract backgroud using top hat in 3D
	 *
	 */

	id = getImageID();
	run("32-bit");
	run("Gaussian Blur 3D...", "x=1 y=1 z=.5");

	run("Duplicate...", "duplicate");
	bg = getImageID();
	run("Minimum 3D...", "x=3 y=3 z=2");
	run("Maximum 3D...", "x=3 y=3 z=2");

	imageCalculator("subtract 32-bit stack", id, bg);
	selectImage(bg); close();

	return id;
}

function findZSlice(channel) {
	/*
	 * Find the brightest slice index
	 *
	 */
	 
	Stack.getDimensions(width, height, channels, slices, frames);
	
	values = newArray(slices);
	max = 0;
	the_slice = 0;
	for (slice = 1; slice <= slices; slice++) {
		Stack.setPosition(channel, slice, 1);
		val = getValue("Mean");
		if ( val > max) {
			max = val;
			the_slice = slice;
		}
	}
	return the_slice;
}

function monoPanel(fig, img, zslice, panel, channel, min, max, lut) {
	selectImage(img);
	// extract middle slice
	Stack.setPosition(channel, zslice, 1);
	// adjust contrast
	setMinAndMax(min, max);
	// set the color mode
	Stack.setDisplayMode("color");
	run(lut);
	// copy
	run("Copy");
	// paste
	selectImage(fig);
	makeRectangle(panel*width, 0, width, height);
	run("Paste");
}

function compositePanel(fig, img, zslice, panel, ch1, min1, max1, lut1, ch2, min2, max2, lut2) {
	selectImage(img);
	// get a composite image with 2 channels
	Stack.setChannel(ch1);
	setMinAndMax(min1, max1);
	run(lut1);
	Stack.setChannel(ch2);
	setMinAndMax(min2, max2);
	run(lut2);
	Stack.setPosition(ch1, zslice, 1);
	Stack.setDisplayMode("composite");
	setActiveChannel(ch1, ch2);
	run("Duplicate...", " ");
	tmp1 = getImageID();
	// convert it to rgb
	run("RGB Color");
	tmp2 = getImageID();
	// copy
	run("Copy");
	// paste in the figure
	selectImage(fig);
	makeRectangle(panel*width, 0, width, height);
	run("Paste");
	// clean up
	selectImage(tmp1); close();
	selectImage(tmp2); close();
}

function histogramPanel(fig, img, ch1, min1, max1, ch2, min2, max2) {
	selectImage(img);
	Stack.getDimensions(width, height, channels, slices, frames);
	h = hist2d(img, ch1, min1, max1, ch2, min2, max2);
	selectImage(h);
	// rescale the histogram to match other panels
	run("Size...", "width="+width+" height="+height+" depth=1 constrain average interpolation=Bilinear");
	// copy
	run("Copy");
	// paste
	selectImage(fig);
	makeRectangle(3*width, 0, width, height);
	run("Paste");
	// clean up
	selectImage(h); close();
}

function printStats(img, channel) {
	selectImage(img);
	run("Duplicate...", "duplicate channels="+channel);
	Stack.getStatistics(voxelCount, mean, min, max, stdDev);
	print("Channel", channel, "range:[",min,"-",max,"]");
}

function makeFigure(img, ch1, min1, max1, lut1, ch2, min2, max2, lut2) {
	/*
	 * Make a colocaliation figure
	 *
	 * Displays the brighest plane of the ch2
	 *
	 * img (int): window id
	 * ch1 (int): channel 1 index
	 * min1 (float)
	 * ch2 (int): channel 2 index
	 * 
	 *
	 */

	 Stack.getDimensions(width, height, channels, slices, frames);
	
	 // select z slice
	 zslice = findZSlice(ch2);
	 print("slice index", zslice);

	 // create the final figure canvas with 4 panels
	 newImage(getTitle, "RGB white", 4*width, height, 1);
	 fig = getImageID();

	 selectImage(img);
	 

	 // for each select channels add a mono panel
	 monoPanel(fig, img, zslice, 0, ch1, min1, max1, lut1);
	 monoPanel(fig, img, zslice, 1, ch2, min2, max2, lut2);
	 
	 // add the composite image
	 compositePanel(fig, img, zslice, 2, ch1, min1, max1, lut1, ch2, min2, max2, lut2);

	 // add a 2D histogram
	 histogramPanel(fig, img, ch1, min1, max1, ch2, min2, max2);
	 selectImage(fig);
	 run("Select None");
	 return fig;
}


start_time = getTime();
run("Close All");
setBatchMode("hide");

// open the image
print(path);
open(path);

// get the ID of the window
img = getImageID();

//preprocess(img);

printStats(img, ch1);
printStats(img, ch2);

// create the colocalization figure
makeFigure(img, ch1, min1, max1, lut1, ch2, min2, max2, lut2);

// save the figure along side the original file
if (savepng) {
	opath = File.getDirectory(path) + File.separator + File.getNameWithoutExtension(path) + ".png";
	print("Saving PNG\n",opath);
	saveAs("PNG", opath);
} else {
	setBatchMode("show");
}
end_time = getTime();
print("Done in " + (end_time - start_time)/1000 + " seconds");
