/*
 * GUV_Background.ijm
 * 
 * Measure average mean intensity over all ROIs for each channel and time frame.
 * 
 * Jerome Boulanger
 */

Stack.getDimensions(width, height, channels, slices, frames);

for (frame = 1; frame <= frames; frame++) {
	for (channel = 1; channel <= channels; channel++) {
		val = 0;
		for (i = 0; i < roiManager("count"); i++) {
			roiManager("select",i);
			Stack.setPosition(channel, 1, frame);
			val += getValue("Mean"); 
		}
		val /= roiManager("count");
		setResult("BG-ch"+(channel), frame-1, val);
	}
}
