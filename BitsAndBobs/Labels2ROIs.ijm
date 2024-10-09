//@String(label="mode", choices={"label","channel"}) colormode
//@Boolean(label="overlay", value=False) asoverlay
/*
 * Convert an image of labels to ROIs
 * 
 * Usage
 * -----
 * Open a label
 * Open the macro and press run
 * 
 * ROI are names as ROI-label-channel-frame and colored by label
 * 
 * Jerome Boulanger 2024
 */
function createLabels() {
	/*Create a test image*/
	width = 320;
	height = 240;
	frames = 10;
	channels = 2;
	labels = 4;
	newImage("HyperStack", "8-bit grayscale-mode", width, height, channels, 1, frames);
	print(labels);
	for (label = 1; label <= labels; label++) {					
		for (channel = 1; channel <= channels; channel++) {						
			x = label/(labels+1) * width;
			y = channel/(channels+1) * height;
			for (frame = 1; frame <= frames; frame++) {
				x += 2 * random - 1;
				y += 2 * random - 1;
				//print(label, channel, frame,x,y);
				Stack.setPosition(channel, 1, frame);				
				setColor(label);
				fillOval(x,y,20,20);		
			}
		}
	}
	for (channel = 1; channel <= channels; channel++) {	
		Stack.setChannel(channel);
		setMinAndMax(0, labels);		
		run("glasbey on dark");
	}
}

function labels2roi(colormode) {
	/*
	 * Convert the current label image to ROI
	 * 
	 * Args:
	 *  colormode: "label" or "channel"	
	 * 
	 */
	setBatchMode("hide");
	iscomposite = false;
	if (nSlices > 1) {
		Stack.getDisplayMode(mode);
		if (mode == "Composite") {
			iscomposite = true;
			Stack.setDisplayMode("grayscale");
		}
	} 
	getDimensions(width, height, channels, slices, frames);
	colors = newArray("red","green","blue","yellow","magenta","orange");
	for (frame = 1; frame <= frames; frame++) {
		for (channel = 1; channel <= channels; channel++) {
			if (nSlices > 1) {
				Stack.setPosition(channel, 0, frame);
			}
			labels = getValue("Max")+1;
			for (label = 1; label <= labels; label++) {
				setThreshold(label-0.1, label+0.1);
				run("Create Selection");
				if(Roi.size>0){					
					Roi.setName("ROI-"+label+"-"+channel+"-"+frame);
					Roi.setPosition(0, 0, frame);
					if (colormode == "label") {
						Roi.setStrokeColor(colors[(label-1)%colors.length]);
					} else {
						Roi.setStrokeColor(colors[(channel-1)%colors.length]);
					}									
					roiManager("add");					
				}
			}
		}
	}
	resetThreshold;
	run("From ROI Manager");
	if (iscomposite) {
		Stack.setDisplayMode("Composite");
	}
	setBatchMode("exit and display");
}


// create a test image if no image are opened
if (nImages==0) createLabels();

// convert the label image to ROI
labels2roi(colormode);

// convert the roi to overlays
if (!asoverlay) { 
	run("To ROI Manager");
}