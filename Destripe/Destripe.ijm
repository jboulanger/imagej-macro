#@String(label="Orientation", choices={"Horizontal","Vertical"}) orientation
#@Integer (label="Bandwidth", value=5) bandwidth
#@Integer (label="Minimum freqency", value=30) min_frequency

function createTestImage() {
	/*
	 * Create a test image with stripes
	 */
	run("Confocal Series");
	run("Macro...", "code=v=v*(0.5+0.1*cos(2*3.14159*x/5))");		
}

function destripe_image(plane, orientation, band,  freq) {
	/*
	 * Destripe an 2D image in Fourier space
	 * 
	 * Args:
	 * 	orientation: orientation of the stripes
	 * 	band: width of the band 
	 * 	freq: low frequency to keep
	 * 
	 */
	run("FFT");
	fft = getImageID();
	getDimensions(width, height, channels, slices, frames);
	if (orientation=="Horizontal") {
		makeRectangle((width - band)/2, 0, band, height);
		run("Cut");
		makeRectangle((width - band)/2, 0, band, height);
		run("Cut");
	} else {
		makeRectangle(0, (height - band)/2, width/2-freq, band);
		run("Cut");
		makeRectangle(width/2+freq, (height - band)/2, width/2-freq, band);
		run("Cut");
	}
	run("Cut");
	run("Inverse FFT");
	ifft = getImageID();
	selectImage(fft);
	close();
	selectImage(plane); 
	close();
	return ifft;
}


function destripe_hyperstack(orientation, band, freq) {
	/*
	 * Destripe the hyperstack
	 * 
	 * Args:
	 * 	orientation: orientation of the stripes
	 * 	band: width of the band 
	 * 	freq: low frequency to keep
	 */
	img = getImageID();
	Stack.getDimensions(width, height, channels, slices, frames);
	for (frame = 1; frame <= frames; frame++) {
		for (slice = 1; slice <= slices; slice++) {
			for (channel = 1; channel <= channels; channel++) {
				selectImage(img);
				Stack.setPosition(channel, slice, frame);
				run("Duplicate...", "duplicate channels="+channel+" slices="+slice+" frame="+frame);
				plane = getImageID();
				result = destripe_image(plane, orientation, band,  freq);
				selectImage(result);
				run("Copy");
				selectImage(img);
				run("Paste");
				selectImage(result);
				close();
			}
		}
	}
}


function appendToFilename(filename, str) {
// append str before the extension of filename
	i = lastIndexOf(filename, '.');
	if (i > 0) {
		name = substring(filename, 0, i);
		ext = substring(filename, i, lengthOf(filename));
		return name + str + ext;
	} else {
		return filename + str
	}
}

function main() {
	// main entry point
	
	// If there is no image create a test image
	if (nImages == 0) {
		createTestImage();
		orientation = "Vertical";
		bandwidth = 15;
		min_frequency = 20;
		
	}
	setBatchMode("hide");
	run("Duplicate...", "title="+appendToFilename(getTitle,"-destriped")+" duplicate");
	destripe_hyperstack(orientation, bandwidth, min_frequency);
	setBatchMode("exit and display");
}

main();