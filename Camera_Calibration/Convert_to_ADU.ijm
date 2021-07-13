/*
 * Open a calibration file and select the corresponding system
 * 
 * The calibration file is created by running "Analyse Camera Noise"
 * 
 * Jerome Boulanger 2018
 */
macro "convert to ADU" {
	systems = newArray(nResults);
	for ( i = 0; i < systems.length; i++) {
		systems[i] = "" + i + " - " + getResultString("System",i);
	}
	Dialog.create("convert to ADU");
	Dialog.addChoice("System",systems);
	Dialog.addNumber("Binning",1);
	Dialog.addNumber("Display max", 20);
	Dialog.show();
	sys = Dialog.getChoice();
	binning = Dialog.getNumber();
	max_display = Dialog.getNumber();	
	row =  -1;
	for (i = 0; i < systems.length && row == -1; i++) {
		if (matches(systems[i], sys)) {
			row = i;
		}
	}
	gain = getResult("A/D Gain", row);
	offset = getResult("DC offset", row);
	edc = getResult("Noise offset (calc.)", row);
	Qe = getResult("Qe", row);	
	if (nSlices>1) {
		Stack.getStatistics(voxelCount, mean, min, max, stdDev);
	} else {
		getStatistics(area, mean, min, max, stdDev, histogram);
	}
	
	psnr = 20*log(tga(max,gain,edc)-tga(offset,gain,edc))/log(10.0);	
	title = getTitle;
	
	run("Duplicate...", "title=["+appendToFilename(title,"-photons")+"] duplicate");
    run("Subtract...", "value="+offset+" stack");
	run("Divide...", "value="+gain+" stack");
	
	if (nSlices>1) {
		Stack.getStatistics(voxelCount, mean, min, max, stdDev);
	} else {
		getStatistics(area, mean, min, max, stdDev, histogram);
	}
	setMinAndMax(0, max_display);
	print("ADU convertion for " + title);
	print("  gain:"+gain);
	print("  offset:"+offset);
	print("  max photon:"+max);
	print("  avg photon:"+mean);
	print("  psnr:"+psnr+"dB");
}

function tga(z,g,edc) {
	return 2/g * sqrt(maxOf(0,g*z+3/8*g*g+edc));
}

// append str before the extension of filename
function appendToFilename(filename, str) {
	i = lastIndexOf(filename, '.');
	if (i > 0) {
		name = substring(filename, 0, i);
		ext = substring(filename, i, lengthOf(filename));
		return name + str + ext;
	} else {
		return filename + str
	}
}