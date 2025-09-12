#@File(label="Input file") input_path
#@Integer(label="Mitochondria channel",value=1) channel_mito
#@Integer(label="ER channel",value=2) channel_er
#@File(label="Mito classifier") class_mito
#@File(label="ER classifier") class_er
#@Boolean(label="Save masks",value=True) save_masks

/*
 * Segment and compute overlap of masks and distance between masks
 * 
 * Open a single z stack
 */


function subtract_background(id, radius_um) {
	selectImage(id);
	getVoxelSize(dx, dy, dz, unit);
	x = 5;maxOf(1,radius_um/dx);
	y = 5;maxOf(1,radius_um/dy);
	z = 2;maxOf(1,radius_um/dz);
	print("background shape ",x,y,z);
	run("Duplicate...", "title=bg duplicate");
	bg = getImageID();
	run("Gaussian Blur 3D...", "x="+x+" y="+y+" z="+z);
	run("Minimum 3D...", "x="+x+" y="+y+" z="+z);
	run("Maximum 3D...", "x="+x+" y="+y+" z="+z);
	imageCalculator("Subtract stack", id, bg);
	selectImage(bg); close();
	selectImage(id);
}

function segment_channel(id, channel) {
	/*
	 * Segment the two channels
	 * 
	 * Parameters
	 *  id: image id
	 *  channel: channel to process
	 *  
	 * Returns
	 *  image id
	 */
	selectImage(id);
	print("Processing channel ", channel);
	background_radius_um = 1;
	getVoxelSize(dx, dy, dz, unit);
	run("Duplicate...", "title=C"+channel+" duplicate channels=" + channel);
	run("32-bit");
	run("Square Root", "stack");
	run("Median...", "radius=1 stack");
	subtract_background(getImageID(), background_radius_um)
	// run("Subtract Background...", "rolling=" + background_radius_um / dx + " stack");
	Stack.getStatistics(voxelCount, mean, min, max, stdDev);
	setThreshold(mean + 4 * stdDev, max);
	run("Make Binary", "background=Dark");
	return getImageID();
}

function segment_labkit(id, channel, classifier_path) {
	print("Segmenting channel " + channel + " with "+ classifier_path);
	selectImage(id);
	run("Duplicate...", "title=C"+channel+" duplicate channels=" + channel);
	run("Segment Image With Labkit", "input=[C"+channel+"] segmenter_file="+classifier_path+" use_gpu=false");
	print(getTitle);
	wait(100);
	selectWindow("segmentation of C"+channel);
	id1 = getImageID();
	run("Duplicate...", "duplicate");
	run("Median...", "radius=1");
	setThreshold(1, 1);
	run("Make Binary", "background=Dark");
	id2 = getImageID();
	selectImage(id1); close();
	return id2;
}

function getVolumeMask(id) {
	selectImage(id);
	Stack.getStatistics(voxelCount, mean, min, max, stdDev);
	return voxelCount * mean / 255;
}


setBatchMode("hide");
run("Close All");

run("Bio-Formats Importer", "open=["+input_path+"] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
run("Select None");
original = getImageID();

// mask_mito = segment_channel(original, channel_mito);
mask_mito = segment_labkit(original, channel_mito, class_mito);
title_mito = getTitle();
run("Distance Transform 3D");
distance_mito = getImageID();

// mask_er = segment_channel(original, channel_er);
mask_er = segment_labkit(original, channel_er, class_er);
title_er = getTitle();
run("Distance Transform 3D");
distance_er = getImageID();

imageCalculator("And create stack", mask_mito, mask_er);
mask_mito_er = getImageID();
title_mito_er = getTitle();

// compute the distance from mito to er
imageCalculator("Multiply stack", distance_er, mask_mito);
selectImage(distance_er);
Stack.getStatistics(voxelCount, mean_distance_mito_er, min, max, stdDev);
mean_distance_mito_er /= 255;

// compute the distance from er to mito
imageCalculator("Multiply stack", distance_mito, mask_er);
selectImage(distance_mito);
Stack.getStatistics(voxelCount, mean_distance_er_mito, min, max, stdDev);
mean_distance_er_mito /= 255;


row = nResults;
volume_mito = getVolumeMask(mask_mito);
volume_er = getVolumeMask(mask_er);
volume_mito_er = getVolumeMask(mask_mito_er);


setResult("Name", row, File.getName(input_path));
setResult("Volume Mito", row, volume_mito);
setResult("Volume ER", row, volume_er);
setResult("Volume Mito x ER", row, volume_mito_er);
setResult("Volume ratio [Mito x ER]:Mito", row, volume_mito_er / volume_mito);
setResult("Volume ratio [Mito x ER]:ER", row, volume_mito_er / volume_er);
setResult("Mean distance Mito:ER", row, mean_distance_mito_er);
setResult("Mean distance ER:Mito", row, mean_distance_er_mito);

run("Merge Channels...", "c1=["+title_mito+"] c2=["+title_er+"] c3=["+title_mito_er+"] create");

output_path = File.getParent(input_path) + File.separator + File.getNameWithoutExtension(input_path) + "-masks.tif";
print("Saving mask as \n"+output_path);
saveAs("TIFF", output_path)

run("Z Project...", "projection=[Max Intensity]");
run("Scale Bar...", "width=2 height=2 horizontal bold overlay");
output_path = File.getParent(input_path) + File.separator + File.getNameWithoutExtension(input_path) + "-masks.png";
print("Saving mip as \n"+output_path);
saveAs("PNG", output_path);

selectImage(original);
run("Z Project...", "projection=[Max Intensity]");
Stack.setDisplayMode("composite");
Stack.setChannel(channel_mito);
run("Enhance Contrast", "saturated=0.35");
Stack.setChannel(channel_er);
run("Enhance Contrast", "saturated=0.35");
run("Scale Bar...", "width=2 height=2 horizontal bold overlay");
output_path = File.getParent(input_path) + File.separator + File.getNameWithoutExtension(input_path) + "-view.png";
print("Saving view as \n"+output_path);
saveAs("PNG", output_path);


setBatchMode("exit and display");
