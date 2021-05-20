#@ Boolean(label="Correct",value=false) correct
/*
 * GUV_Background.ijm
 * 
 * Measure average mean intensity over all ROIs for each channel and time frame.
 * 
 * correct: add a corrected column or estimate the background
 * 
 * Jerome Boulanger 2021
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
		if (correct) {
			header = split(Table.headings,"\t");
			for (i = 0; i < header.length; i++){
				if (matches(header[i], ".*-ch"+channel)) {
					for (row = 0; row < nResults; row++) {
						I = getResult(header[i], row);
						setResult(header[i]+"-corrected", row, I-val);
					}
				}
			}
		} else {
			setResult("BG-ch"+(channel), frame-1, val);
		}
	}
}
