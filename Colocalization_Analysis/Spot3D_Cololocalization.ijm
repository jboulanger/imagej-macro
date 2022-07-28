#@File (label="Input file") ipath
#@String(label="Channel indices", values="1,2,3") channel_indices_str
#@String(label="Channel names", values="Bassoon,GluA2,PSD95") channel_names_str
#@Float(label="Blob size [um]", value=0.125) blob_size_um
#@Float(label="Specificity", value=8) blob_pfa
#@Boolean(label="Save and Close",value=true) dosave
#@File (label="Output directory", style="directory") ofolder

/*
 * 3D Spots colocalization
 * 
 * - Detect spots in the selected channels
 * - Save the spot centroids in a table
 * 
 * 
 * Jerome Boulanger for Viktor & Imogen 2022
 */


run("Close All");

setBatchMode("hide");

start_time = getTime();

channel_indices = parseCSVInt(channel_indices_str);
channel_names = parseCSVString(channel_names_str);
name = File.getName(ipath);

tbl = "coordinates.csv";
if (isOpen(tbl)){selectWindow(tbl);run("Close");}
Table.create(tbl);

id = loadImage(ipath);

labels = segmentImage(id, channel_indices, blob_size_um, blob_pfa);

appendCentroidToTable(tbl, name, channel_names, channel_indices, labels);

computeDistances(tbl, channel_names, channel_indices);

merge = mergeLabels(labels, channel_names, name);

if (dosave) saveOutput(ofolder, tbl, merge, ipath);

end_time = getTime();

print("Elapsed Time " + (end_time - start_time)/ 1000 + " seconds.");

setBatchMode("exit and display");

function loadImage(path) {
	// load the image and return its id
	print("Opening image");
	print(path);
	run("Bio-Formats Importer", "open=["+path+"] color_mode=Composite rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
	return getImageID();	
}

function saveOutput(folder, tbl, merge, path) {
	// save the table and label merged image in the folder
	selectWindow(tbl);
	opath = folder + File.separator + File.getNameWithoutExtension(path) + ".csv";
	print("Saving table to file");
	print(opath);
	Table.save(opath);
	// save the color merged label
	selectImage(merge);	
	opath = folder + File.separator + File.getNameWithoutExtension(path) + "-label.tif";
	print("Saving labels to file");
	print(opath);
	saveAs("TIFF", opath);	
}

function segmentImage(id, channels, size, pfa) {	
	// segment each channels from image id
	labels = newArray(channels.length);
	for (k = 0; k < channels.length; k++) {		
		labels[k] = segment3DBlobs2(id, channels[k], size, pfa, 0);
	}	
	return labels;
}

function mergeLabels(labels, channel_names, name) {
	// merge labels and return the id of the image
	str = "";
	for (k = 0; k< labels.length; k++) {
		selectImage(labels[k]);
		rename(channel_names[k]);
		str += "c"+(k+1)+"=["+channel_names[k]+"] ";
	}
	run("Merge Channels...", str + "create");
	rename(File.getNameWithoutExtension(name)+"-label.tif");
	
	for (k = 0; k< labels.length; k++) {
		Stack.setChannel(k+1);
		getLut(r, g, b);		
		setColor(r[128],g[128],b[128]);
		Overlay.drawString(channel_names[k], 10, 10+12*k);
		Overlay.setPosition(0);
		Overlay.show();
		setMinAndMax(0, 1);
	}
	return getImageID();
}

function appendCentroidToTable(tbl, name, channel_names, channel_indices, labels) {
	// append centroids to the table for each label image	
	for (k = 0; k < labels.length; k++) {
		selectImage(labels[k]);		
		run("Analyze Regions 3D", "volume centroid surface_area_method=[Crofton (13 dirs.)] euler_connectivity=6");
		wait(500);		
		volume = Table.getColumn("Volume");
		id = Table.getColumn("Label");
		x = Table.getColumn("Centroid.X");	
		y = Table.getColumn("Centroid.Y");	
		z = Table.getColumn("Centroid.Z");
		run("Close");
		selectWindow(tbl);
		row0 = Table.size;
		for (i = 0; i < id.length; i++) {	
			row = row0 + i;
			Table.set("Image Name", row, name);
			Table.set("Index", row, row);			
			Table.set("Channel Name", row, channel_names[k]);			
			Table.set("Channel Index", row, channel_indices[k]);	
			Table.set("Label", row, id[i]);
			Table.set("Volume", row, volume[i]);
			Table.set("X", row, x[i]);
			Table.set("Y", row, y[i]);
			Table.set("Z", row, z[i]);	
			Table.set("C", row, k);	
		}
	}
}

function computeDistances(tbl, channel_names, channel_indices) {
	// compute distances between centroid spots 
	selectWindow(tbl);
	N = Table.size;
	print(" Number of objects : " + N);
	count = 0;
	// for each spot
	for (i = 0; i < N; i++) {
		count++;
		showProgress(i, N);
		showStatus("computing distances " + i + "/" + N + "~"+count);
		// get the position and channel		
		xi = Table.get("X", i);
		yi = Table.get("Y", i);
		zi = Table.get("Z", i);
		ci = Table.get("C", i);
		// search for the nearest neightbor in each channel
		idx = newArray(channel_indices.length);
		dmin = newArray(channel_indices.length);		
		for (j = 0; j < channel_indices.length; j++) {
			dmin[j] = 1e6;
			idx[j] = -1;
		}		
		for (j = 0; j < N; j++) {
			cj = Table.get("C",j);				
			if (ci != cj) {							
				dx = Table.get("X", j) - xi;
				dy = Table.get("Y", j) - yi;
				dz = Table.get("Z", j) - zi;		
				dij = sqrt(dx*dx+dy*dy+dz*dz);
				if (dij < dmin[cj]) {				
					dmin[cj] = dij;
					idx[cj] = j;
				}
			}
		}
		// store the NN information in the table
		for (j = 0; j < channel_indices.length; j++) {	
			if (idx[j] > 0) {
				Table.set(channel_names[j] + " NN index", i, idx[j]);
				Table.set(channel_names[j] + " NN distance", i, dmin[j]);			
			} else {
				Table.set(channel_names[j] + " NN index", i, -1);
				Table.set(channel_names[j] + " NN distance", i, NaN);			
			}
		}
	}
	Table.update;	
}

function parseCSVString(csv) {
	str = split(csv, ",");
	for (i = 0; i < str.length; i++) { 
		str[i] =  String.trim(str[i]);
	}
	return str;
}

function parseCSVInt(csv) {
	a =  split(csv, ",");
	for (i = 0; i < a.length; i++){ 
		a[i] =  parseInt(String.trim(a[i]));
	}
	return a;
}


function segment3DBlobs(id, channel, size, pfa, maskid) {
	selectImage(id);name=getTitle();	
	print("Segment3DBlob "+getTitle()+" channel:"+channel+", size:"+size+"um, pfa:"+pfa);
	getVoxelSize(dx, dy, dz, unit);	
	run("Select None");
	sigma1 = size;
	sigma2 = 2 * size;	
			
	run("Duplicate...", "duplicate channels="+channel);	
	id1 = getImageID();rename("id1");
	run("32-bit");
	run("Square Root","stack");
		
	run("Gaussian Blur 3D...", "x="+sigma1/dx+" y="+sigma1/dy+" z="+sigma1/dz);
	run("Duplicate...", "duplicate");
	id2 = getImageID();
	run("Gaussian Blur 3D...", "x="+sigma2/dx+" y="+sigma2/dy+" z="+sigma2/dz);	
	imageCalculator("Subtract 32-bit stack",id1,id2);
	selectImage(id2);close();		

	// threshold score map
	selectImage(id1);
	threshold(pfa);	
	run("Make Binary", "method=Default background=Dark black");	
	rename("mask");
	
	run("Connected Components Labeling", "connectivity=26 type=float");
	id3 = getImageID();
	//run("Label Size Filtering", "operation=Lower_Than size=1000");
	//id4  =getImageID();
	
	selectImage(id1);close();
	//selectImage(id3);close();
	selectImage(id3); rename("Spots in "+name);
	return id3;	
}

function segment3DBlobs2(id, channel, size, pfa, maskid) {
	selectImage(id);name=getTitle();	
	print("Segment3DBlob "+getTitle()+" channel:"+channel+", size:"+size+"um, pfa:"+pfa);
	getVoxelSize(dx, dy, dz, unit);	
	run("Select None");
				
	// difference of Gaussian
	sigma1 = size;
	sigma2 = 2 * size;	
	run("Duplicate...", "duplicate channels="+channel);	
	id1 = getImageID(); rename("id1");
	run("32-bit");
	run("Square Root","stack");	
	run("Gaussian Blur 3D...", "x="+sigma1/dx+" y="+sigma1/dy+" z="+sigma1/dz);
	run("Duplicate...", "duplicate");
	id2 = getImageID();
	run("Gaussian Blur 3D...", "x="+sigma2/dx+" y="+sigma2/dy+" z="+sigma2/dz);	
	imageCalculator("Subtract 32-bit stack",id1,id2);
	selectImage(id2);close();		

	// local maxima
	selectImage(id1);
	run("Duplicate...", "duplicate");
	// smooth the image first
	id3 = getImageID();
	run("Gaussian Blur 3D...", "x="+sigma2/dx+" y="+sigma2/dy+" z="+sigma2/dz);
	// compute regional max
	run("Regional Min & Max 3D", "operation=[Regional Maxima] connectivity=26");	
	id4 = getImageID();
	// close the smoothed image
	selectImage(id3);close();
	
	// label the  local maximas
	selectImage(id4);
	run("Connected Components Labeling", "connectivity=26 type=float");
	id5 = getImageID();
	rename("label");
	
	// compute the distance from the local maximas image
	selectImage(id4);
	run("Macro...", "code=v=255-v stack");	
	run("Chamfer Distance Map 3D", "distances=[Quasi-Euclidean (1,1.41,1.73)] output=[32 bits]");
	id6 = getImageID();
	rename("distance");
	
	// threshold score map
	selectImage(id1);	
	threshold(pfa);
	run("Make Binary", "method=Default background=Dark black");	
	rename("mask");
	
	// separate touching objects
	run("Marker-controlled Watershed", "input=distance marker=label mask=mask compactness=0 binary calculate use");
	id7 = getImageID();
		// filter by size
	v0 = round(100 * 4/3*PI*size*size*size/(dx*dy*dz));	
	run("Label Size Filtering", "operation=Lower_Than size="+v0);
	id8 = getImageID();
	run("Label Size Filtering", "operation=Greater_Than size="+7);
	id9 = getImageID();
	
	selectImage(id1);close();
	selectImage(id4);close();
	selectImage(id5);close();
	selectImage(id6);close();
	selectImage(id7);close();
	selectImage(id8);close();
	selectImage(id9); rename("Spots Ch"+channel+" in "+name);
	return id9;	
}

function threshold(logpfa) {	
	 // Define a threshold based on a probability of false alarm 
	 st = laplaceparam();
	 pfa = pow(10,-logpfa);	 
	 if (pfa > 0.5) {
	 	t = st[0] + st[1]*log(2*pfa);
	 } else {
	 	t = st[0] - st[1]*log(2*pfa);
	 }	 
	 setThreshold(t, st[2]);
}

function laplaceparam() {
	/** Compute laplace distribution parameters 
	 *  Return center, scale and maximum value
	**/
	nbins = 500;
	values = getStackHistogramValues(nbins);
	counts =  getStackHistogramCounts(values);
	Stack.getDimensions(width, height, channels, slices, frames);	
	n = width*height*slices*channels*frames;
	c = 0;
	for (i = 0; i < counts.length && c <= n/2; i++) {
		c = c + counts[i];
	}
	m = values[i-1];
	b = 0;
	for (i = 0; i < counts.length; i++){
		b = b + abs(values[i] - m) * counts[i];
	}
	max =  2 * values[values.length-1] - values[values.length-2];
	return newArray(m,b/n,max);
}

function getStackHistogramValues(nbins) {	
	Stack.getStatistics(voxelCount, mean, min, max, stdDev);
	values = newArray(nbins);
	for (i = 0; i < nbins; i++) {
		values[i] = min + (max - min) * i / nbins;
	}
	return values;
}

function getStackHistogramCounts(values) {
	Stack.getDimensions(width, height, channels, slices, frames);	
	min = values[0];
	max =  2 * values[values.length-1] - values[values.length-2];
	nbins= values.length;
	hist = newArray(nbins);
	for (slice = 1; slice <= nSlices; slice++) {
		setSlice(slice);
		getHistogram(vals, counts, nbins, min, max);
		for (i=0;i<nbins;i++) {
			hist[i] = hist[i] + counts[i];
		}
	}
	return hist;
}
