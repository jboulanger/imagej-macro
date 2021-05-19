/*
 * Analyse the intensity of GUV over time
 * 
 * Data are 2D image sequence with two channels.
 * The spherical GUVs apprear as bright rings in the images and 
 * their motion is very limited.
 * 
 * The signal is segmented using a difference of Gaussian in both channels
 * an OR operation provide a unique mask from the two channels.
 * 
 * The ROI from the user is then utilized to solve the tracking problem 
 * and identify each GUV for each other over time
 * 
 * Intensities overtime is then measured for each channel and object and 
 *  save in the result table.
 * 
 * Usage:
 * - Open the original files using Bioformat and select open as hyperstack levaving evering thing else unchecked.
 * - Compute a maximum intensity projection to see where the GUV is located overall frame
 * - Select several ROI (draw region and press t) enclosing each GUV during the all movie
 * - Launch the macro (Run in the editor)
 * 
 * Jerome Boulanger 2016 for Lauren McGinney and Yohei Ohashi
 * Modified in mar 2018 for error of image focus + outlines
 */
 
macro "GUV Intensity" {	
	Stack.getDimensions(width, height, channels, slices, frames);
	Dialog.create("Parameters");
	Dialog.addSlider("Smoothing [%]",0.0,100.0,50.0);
	Dialog.addSlider("Sensitivity [%]",0.0,100.0,50.0);	
	Dialog.addSlider("Width [px]",0.0,20.0,10.0);
	Dialog.addChoice("Measure", newArray("Mean","RawIntDen","Median","IntDen"));
	Dialog.addSlider("First frame",1,frames,1);
	Dialog.addSlider("Last frame",1,frames,frames);	
	Dialog.show();
	smoothing = Dialog.getNumber();
	sensitivity = Dialog.getNumber();	
 	bandwidth = Dialog.getNumber();
 	measure = Dialog.getChoice();
 	first_frame = Dialog.getNumber();
 	last_frame = Dialog.getNumber();	
 	print("GUV intensity ["+smoothing+","+sensitivity+","+bandwidth+","+measure+"]"); 	
	guvIntensity(smoothing, sensitivity, bandwidth, measure, first_frame, last_frame);
	//segmentAllGUVs(smoothing/100*3, 4*(100-sensitivity)/100, bandwidth);
	showGraphs(channels);
}

function guvIntensity(smoothing, sensitivity, bandwidth, measure, first_frame, last_frame) {
 	
	// Check the number of ROI in the ROI manager
	nroi = roiManager("count");
	while (nroi == 0) {
		waitForUser("Please select ROI and add them to the ROI manager");
		nroi = roiManager("count");
	}
	print("Analysis of " + nroi + " regions.");
			
	// Clean up initialization
	run("Set Measurements...", "None redirect=None decimal=3");	
	run("Select None");
	run("Clear Results"); 
	Overlay.remove();
	setBatchMode(true);
	id = getImageID;
	segmentAllGUVs(smoothing/100.0*3.0, 4.0*(100-sensitivity)/100.0, bandwidth, first_frame, last_frame);
	selectImage(id);
	measureIntensities(bandwidth, measure, first_frame, last_frame);
	
	// Clean up created ROIs	
	sel = newArray(roiManager("count") - nroi);	
	for (i = 0; i< sel.length; i++) {
		sel[i] = nroi + i;
	}
	roiManager("Select", sel);
	roiManager("Delete");
	run("Select None");
	
	Stack.setFrame(1);
	Stack.setChannel(1);
	setBatchMode(false);
	print("Done");
}

// Segmentation of all images
function segmentAllGUVs(smoothing, threshold, bandwidth, first_frame, last_frame) {	
	
	id = getImageID;	
	run("Duplicate...", "duplicate channels=1");
	id1 = getImageID;
	
	segment(smoothing, threshold);
	selectImage(id);
	run("Duplicate...", "duplicate channels=2");
	id2 = getImageID;	
	segment(smoothing, threshold);
	
	imageCalculator("OR stack", id1, id2);	
	getOutlines(bandwidth);	
	setThreshold(128,255);
		
	for (t = first_frame; t <= last_frame; t++) {
		Stack.setFrame(t);
		run("Create Selection");
		roiManager("Add");		
	}
	selectImage(id1);close();
	selectImage(id2);close();
	selectImage(id);
	//selectWindow("ROI Manager");	
	}
}

function segment(scale, lambda) {	
	id1 = getImageID;
	run("32-bit");
	run("Square Root", "stack");
	run("Gaussian Blur...", "sigma="+scale+" stack");
	run("Duplicate...", "duplicate ");
	id2 = getImageID;
	run("Gaussian Blur...", "sigma="+3*scale+" stack");	
	imageCalculator("Subtract stack 32-bit", id1, id2);
	selectImage(id2); close();
	selectImage(id1);	
	Stack.getStatistics(voxelCount, mean, min, max, std)
	threshold = mean + lambda * std;
	setThreshold(threshold, max);	
	run("Convert to Mask", "method=Default background=Default");			
	smoothContour(20,5,3);	
	}
}

function getOutlines(size) {
	id = getImageID;
	run("Duplicate...", "duplicate ");
	id1 = getImageID;	
	run("Fill Holes", "stack");	
	smoothContour(20,5,3);	
	run("Duplicate...", "duplicate ");
	id2 = getImageID;
	run("Minimum...", "radius="+size+" stack");	
	imageCalculator("XOR stack", id1, id2);	
	selectImage(id2); close();
	imageCalculator("AND stack", id, id1);	
	selectImage(id1); close();
	selectImage(id);	
}

function smoothContour(iter, r1, r2) {
	for (i = 0; i < iter; i++) {
		run("Maximum...", "radius="+r1+" stack");
		run("Minimum...", "radius="+r1+" stack");
		run("Median...", "radius="+r2+" stack");
	}
}

// Measure intensities for each ROI and frame
function measureIntensities(bandsize, measure, first_frame, last_frame)  {
	has_empty_rois = false;
	getDimensions(width, height, channels, slices, frames);
	for (i = 0; i < nroi; i++) {
		roiManager("select", i);
		name = Roi.getName();	
		for (t = first_frame; t <= last_frame; t++) {
			Stack.setFrame(t);
			// create a region in the ROI from the mask			
			roiManager("select", newArray(i, nroi + t - 1));
			roiManager("AND");		
			getSelectionBounds(x, y, w, h);
			if ((w != width) || (h != height)) {
				// measure intensities
				setResult("Frame", t-1, t);
				for (c = 1; c <= channels; c++) {
					Stack.setChannel(c);
					Stack.setFrame(t);					
					List.setMeasurements;
					setResult("ROI"+(i+1) + "-ch"+c, t-1, d2s(List.getValue(measure),2));		
					sn = getSliceNumber();
					// add an overlay					
					run("Add Selection...", "stroke=#9999ff width=1");
					Overlay.setPosition(c, 1, t);
					Overlay.add;		
				}
			} else {	
				has_empty_rois = true;
				// Error: the Roi is empty, we just put -1 in the table
				for (c = 1; c <= channels; c++) {
					setResult(name + "-ch"+c, t-1, -1);
				}
			}
		}
	}
	updateResults();
	if (has_empty_rois)  {
		print("*** Warning ***\nUnable to segment all objects.\nSetting values to -1 in the result table.\nPlease try to modify the sensitivity\n***********\n");
	}
}

function showGraphs(nchannels) {
	// Load the table and create graphs
	nroi = roiManager("count");
	Plot.create("Graph","Frame","Intensity")
	legend="";
	colors = getLutHexCodes("Ice",nchannels*nroi*2);
	Array.print(colors);
	mini=255;
	maxi=0;
	type = newArray("ROI");
	for (i = 1; i <= nroi; i++) {
		for (c = 1; c <= nchannels; c++) {
			for (k = 0; k < type.length; k++) {
				T = Table.getColumn("Frame");
				I = Table.getColumn(""+type[k]+i+"-ch"+c);
				Array.getStatistics(I, min, max);
				mini = minOf(min,mini);
				maxi = maxOf(max,maxi);
				legend += ""+type[k]+i+"-ch"+c+"\t";
				Plot.setColor(colors[k+2*((c-1)+nchannels*(i-1))]);
				Plot.setLineWidth(2);
				Plot.add("line", T, I,""+type[k]+i+"-ch"+c);
			}
		}
	}
	Plot.setLimits(0, T[T.length-1], mini, maxi);
	Plot.setLegend(legend, "options");
}

function getLutHexCodes(name,n) {
	//setBatchMode("hide");
	newImage("stamp", "8-bit White", 10, 10, 1);
   	run(name);
   	getLut(r, g, b);
   	close();
   	r = Array.resample(r,n);
   	g = Array.resample(g,n);
   	b = Array.resample(b,n);
	dst = newArray(r.length);
	for (i = 0; i < r.length; i++) {
		dst[i] = "#"+dec2hex(r[i]) + dec2hex(g[i]) + dec2hex(b[i]);
	}   
	//setBatchMode("show");
	return dst;	
}

function dec2hex(n) {
	// convert a decimal number to its heximal representation
	codes = newArray("0","1","2","3","4","5","6","7","8","9","A","B","C","D","E","F");
	k1 = Math.floor(n/16);
	k2 = n - 16 * Math.floor(n/16);
	return codes[k1]+codes[k2]; 
}

