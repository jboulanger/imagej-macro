#@File (label="Input file (csv)") ipath_csv
#@String(label="Channel names", values="bassoon,glua2,psd95") channel_names_str_order
#@Float(label="Distance threshold [um]",value=1) distance_threshold
#@Boolean(label="Visualization", value=false) clusters_visualization

/*
 * Analyze Spot3D Distances
 * 
 * This script analyses the results of the macro Spot3D_Distances
 * 
 * Find chains of markers in the order indicated by the names such that the distance between
 * each pair is below a threshold.
 * 
 * Output:
 * table overlaps.csv: number of spots and number of spots < distance, mean distance for d<threshold
 * table clusters.csv: report spots chained clusters
 * 
 * The name of the channels need to matche the one of the input table
 * 
 * If processing several files with "batch", results will be appended in each table.
 * 
 */

var file_count; // index of file (global var)

// skip sliently the file if not a csv file
if (!endsWith(ipath_csv, ".csv")) {exit();}

setBatchMode("hide");

start_time = getTime();

channel_names = parseCSVString(channel_names_str_order);

// Open the csv file
print("Loading\n" + ipath_csv);
open(ipath_csv);
tbl1 = Table.title;
file_count++;

print("Chained clusters");
tbl3 = "clusters.csv";
if (clusters_visualization) {
	print("Creating new table " + tbl3);
	Table.create(tbl3); // erase the table content
} else {		
	if (!isOpen(tbl3)) {		
		print("Create new table" + tbl3);
		Table.create(tbl3);  // keep the previous lines
	} else {
		print("Reusing table " + tbl3);
	}
}

indices = identifyChainedClusters(tbl1, tbl3, channel_names, distance_threshold);
nchained = indices.length;

print("Overlap measurements");
tbl2 = "overlaps.csv";
if (clusters_visualization) {		
	print("Creating new table " + tbl2);
	Table.create(tbl2); // erase the table content
} else {
	if (!isOpen(tbl2)) { 		
		print("Creating new table " + tbl2);
		Table.create(tbl2); // keep the previous lines
	} else {
		print("Reusing table " + tbl2);
	}
} 
computeOverlap(tbl1, tbl2, channel_names, distance_threshold, nchained);


selectWindow(tbl1); run("Close");



if (clusters_visualization) {	
	folder = File.getParent(File.getDirectory(ipath_csv));
	displayChainedCluster(folder, tbl3, indices, distance_threshold);
}

setBatchMode("exit and display");

end_time = getTime();

print("Elapsed Time " + (end_time - start_time)/ 1000 + " seconds.\n");



function displayChainedCluster(folder, clusters, idx, distance_threshold) {
	print("Visualization");	
	selectWindow(clusters);
	name = Table.getString("Image Name", idx[0]);	
	print("Loading");
	print(folder + File.separator + name);
	open(folder + File.separator + name);	
	getVoxelSize(dx, dy, dz, unit);
	r = round(distance_threshold/dx);
	colors  = newArray("Red","Green","Blue");
	for (i = 0 ; i<3; i++) {Stack.setChannel(i+1); run(colors[i]); }
	img = getImageID();
	setColor("yellow");
	Overlay.remove();	
	Stack.getDimensions(width, height, channels, slices, frames);
	print("Image size" + width + "x" + height + "x" + slices );
	for (i = 0; i < idx.length; i++) {
		selectWindow(clusters);		
		x = round(Table.get("X", idx[i]) / dx);
		y = round(Table.get("Y", idx[i]) / dy);
		z = round(Table.get("Z", idx[i]) / dz);		
		print("X:"+x+", Y:"+y+", Z:"+z + ", R:"+r);
		selectImage(img);
		Overlay.drawEllipse(x-r, y-r, 2*r, 2*r);
		Overlay.add();
		Overlay.setPosition(0, z, 0);		
	}	
	Overlay.show();
}

function computeOverlap(src, dst, channel_names, threshold, nchained) {	
	// count the number of spots in each channel
	// count the number of time a nearest neighbor is below the threshold distance
	// measure the average distance to nearest neighbor under the threshold
	// report results in a table
	
	selectWindow(src);
	N = Table.size;
	counts = newArray(channel_names.length);
	nnoverlaps = newArray(channel_names.length * channel_names.length);	
	nndistances = newArray(channel_names.length * channel_names.length);	
	for (i = 0; i < N; i++) {
		selectWindow(src);
		chi = Table.getString("Channel Name", i);
		// find the channel name of this line and update count
		k0 = -1; // channel 
		for (k = 0; k < channel_names.length; k++) {			
			if (matches(chi, channel_names[k])) {
				counts[k] = counts[k] + 1;
				k0 = k;				
			}
		}
		if (k0 >= 0) {
			for (k = 0; k < channel_names.length; k++) if (k != k0) {
				d = Table.get(channel_names[k] + " NN distance", i);						
				if (d < threshold && !isNaN(d)) {
					I = k0 + channel_names.length * k;
					nnoverlaps[I] = nnoverlaps[I] + 1;
					nndistances[I] = nndistances[I] + d;				
				}
			}
		}
	}
	
	// report results in a table
	name = Table.getString("Image Name");
	selectWindow(dst);
	row = Table.size;	
	Table.set("Image Name", row, name);
	Table.set("Image Index", row, file_count);
	Table.set("Threshold", row, threshold);
	Table.set("# chained", row, nchained);
	for (k = 0; k < channel_names.length; k++) {	
		Table.set(channel_names[k], row, counts[k]);
	}
	for (k = 0; k < channel_names.length; k++) {
		for (l = 0; l < channel_names.length; l++) if (k!=l) {
			Table.set("#" + channel_names[k]+"-"+channel_names[l], row, nnoverlaps[k+channel_names.length*l]);			
		}		
	}
	for (k = 0; k < channel_names.length; k++) {
		for (l = 0; l < channel_names.length; l++) if (k!=l) {
			I = k + channel_names.length * l;
			Table.set("Ratio ["+channel_names[k]+":"+channel_names[l]+"]", row, nnoverlaps[I] / counts[k]);
		}		
	}
	for (k = 0; k < channel_names.length; k++) {
		for (l = 0; l < channel_names.length; l++) if (k!=l) {
			I = k + channel_names.length * l;
			Table.set("Distance ["+channel_names[k]+":"+channel_names[l]+"]", row, nndistances[I] / nnoverlaps[I]);
		}		
	}
	Table.update;
}


function identifyChainedClusters(src, dst, channel_names, threshold) {
	// follow a chain of markers based on a threshold
	// reports the indices of the markers, distances and volumes
	selectWindow(dst);
	row0 = Table.size;	
	row = row0;
	selectWindow(src);
	N = Table.size;
	name = Table.getString("Image Name");	
	for (i = 0; i < N; i++) {
		selectWindow(src);		
		chi = Table.getString("Channel Name", i);
		vi = Table.getString("Channel Name", i);
		idx = newArray(channel_names.length);
		distances = newArray(channel_names.length-1);
		volumes = newArray(channel_names.length);
		if (matches(chi, channel_names[0])) { // start of the chain
			idx[0] = i;
			volumes[0] = Table.get("Volume", i);
			x = Table.get("X", i);
			y = Table.get("Y", i);
			z = Table.get("Z", i);
			j = i; // pointer along the chain			
			for (k = 1; k < channel_names.length; k++) {
				d = Table.get(channel_names[k] + " NN distance", j);				
				idx[k] = Table.get(channel_names[k] + " NN index", j);
				distances[k-1] = d;
				volumes[k] = Table.get("Volume", j);			
			}
			Array.getStatistics(distances, dmin, dmax);			
			if (dmax < threshold) {				
				selectWindow(dst);
				row = Table.size;
				Table.set("Image Name", row, name);
				Table.set("Image Index", row, file_count);
				for (k = 0; k < channel_names.length; k++) {
					Table.set(channel_names[k], row, idx[k]);					
				}				
				for (k = 0; k < distances.length; k++) {
					Table.set("Distance ["+channel_names[k]+":"+channel_names[k+1]+"]", row, distances[k]);					
				}				
				for (k = 0; k < channel_names.length; k++) {
					Table.set("Volume [" + channel_names[k] + "]", row, volumes[k]);
				}
				Table.set("X", row, x);
				Table.set("Y", row, y);
				Table.set("Z", row, z);				
			}
			Table.update;
		}
	}
	// generate the list of rows that were added
	idx = newArray(row-row0);
	for (i = 0; i < idx.length; i++) {
		idx[i] = row0 + i;
	}	
	return idx;
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