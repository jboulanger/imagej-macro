//@File(style="open") input

/*
 * Extract time signals from ROI in the roi manager
 * 
 * Jerome Boulanger for Alex Fellow 2023
 */
 
 run("Close All");
 if (isOpen("ROI Manager")) {
 	selectWindow("ROI Manager"); 
 	run("Close");
 }

 start_time = getTime();
 setBatchMode("hide");
 run("TIFF Virtual Stack...", "open=["+input+"]");
 roiManager("Open", replace(input,".tif","_ROI.zip"));
 open(replace(input,".tif","_ROI.zip"));
 Stack.getDimensions(width, height, channels, slices, frames);
 N  = roiManager("count");
 for (frame = 1; frame < frames; frame++) {
 	Stack.setFrame(frame);
 	for (i = 0; i < N; i++) {
 		roiManager("select", i);
 		roiname = Roi.getName;
 		setResult(roiname, frame-1, getValue("Mean"));
 	}
 }
updateResults();
setBatchMode("show");
print("Elapsed time " + (getTime() - start_time)/1000 + " seconds.");
