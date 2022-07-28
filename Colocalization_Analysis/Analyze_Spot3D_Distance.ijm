#@File (label="Input file") ipath
#@String(label="Channel names", values="Bassoon,GluA2,PSD95") channel_names_str
#@Float(label="Distance threshold [um]",value=1) distance_threshold
/*
 * Find chains of markers in the order indicated by the names such that the distance between
 * each pair is below a threshold.
 * 
 * The name of the channels need to matche the one of the input table
 * 
 */


run("Close All");

setBatchMode("hide");

start_time = getTime();

channel_names = parseCSVString(channel_names_str);

open(ipath);
tbl1 = Table.title;

tbl2 = "overlaps.csv";
if (!isOpen(tbl2)){Table.create(tbl2);}

statistics(tbl1, tbl2, channel_names, distance_threshold);

tbl3 = "distances.csv";
if (!isOpen(tbl3)){Table.create(tbl3);}
identifyChainedClusters(tbl1, tbl3, channel_names, distance_threshold);

//selectWindow(tbl1); run("Close");

end_time = getTime();

print("Elapsed Time " + (end_time - start_time)/ 1000 + " seconds.");

setBatchMode("exit and display");


function statistics(src, dst, channel_names, threshold) {	
	selectWindow(src);
	N = Table.size;
	counts = newArray(channel_names.length);
	overlaps = newArray(channel_names.length * channel_names.length);	
	distances = newArray(channel_names.length * channel_names.length);	
	for (i = 0; i < N; i++) {
		chi = Table.getString("Channel Name", i);
		k0 = 1;
		for (k = 0; k < channel_names.length; k++) {			
			if (matches(chi, channel_names[k])) {
				counts[k] = counts[k] + 1;
				k0 = k;
				break;
			}
		}
		for (k = 0; k < channel_names.length; k++) {
			d = Table.getString(channel_names[k] + " NN distance", i);			
			if (d < threshold && !isNaN(d)) {
				I = k0+channel_names.length*k;
				overlaps[I] = overlaps[I] + 1;
				distances[I] = distances[I] + d;
				//I = k+channel_names.length*k0;
				//overlaps[I] = overlaps[I] + 1;
				//distances[I] = distances[I] + d;
			}
		}
	}
	
	// report results in a table
	name = Table.getString("Image Name");
	selectWindow(dst);
	row = Table.size;	
	Table.set("Image Name", row, name);
	Table.set("Threshold", row, threshold);
	for (k = 0; k < channel_names.length; k++) {	
		Table.set(channel_names[k], row, counts[k]);
	}
	for (k = 0; k < channel_names.length; k++) {
		for (l = 0; l < channel_names.length; l++) if (k!=l) {
			Table.set(channel_names[k]+"-"+channel_names[l], row, overlaps[k+channel_names.length*l]);
		}		
	}
	for (k = 0; k < channel_names.length; k++) {
		for (l = 0; l < channel_names.length; l++) if (k!=l) {
			Table.set("ratio ["+channel_names[k]+"-"+channel_names[l]+"]", row, overlaps[k+channel_names.length*l]/counts[k]);
		}		
	}
	for (k = 0; k < channel_names.length; k++) {
		for (l = 0; l < channel_names.length; l++) if (k!=l) {
			Table.set("delta ["+channel_names[k]+"-"+channel_names[l]+"]", row, distances[k+channel_names.length*l] / overlaps[k+channel_names.length*l]);
		}		
	}
	Table.update;
}


function identifyChainedClusters(src, dst, channel_names, threshold) {
	// follow a chain of markers based on a threshold
	// reports the indices of the markers, distances and volumes
	selectWindow(src);
	N = Table.size;
	row = 0;
	for (i = 0; i < N; i++) {
		selectWindow(src);		
		chi = Table.getString("Channel Name", i);
		vi = Table.getString("Channel Name", i);
		idx = newArray(channel_names.length);
		distances = newArray(channel_names.length-1);
		volumes = newArray(channel_names.length);
		if (matches(chi, channel_names[0])) { // start of the chain
			idx[0] = i;
			volumes[0] = Table.getString("Volume", i);
			j = i; // pointer along the chain			
			for (k = 1; k < channel_names.length; k++) {
				d = Table.getString(channel_names[k] + " NN distance", j);				
				idx[k] = Table.get(channel_names[k] + " NN index", j);
				distances[k-1] = d;
				volumes[k] = Table.get("Volume", j);			
			}
			Array.getStatistics(distances, dmin, dmax);			
			if (dmax < threshold) {				
				selectWindow(dst);
				row = Table.size;
				for (k = 0; k < channel_names.length; k++) {
					Table.set(channel_names[k], row, idx[k]);					
				}				
				for (k = 0; k < distances.length; k++) {
					Table.set(channel_names[k]+"-"+channel_names[k+1], row, distances[k]);					
				}				
				for (k = 0; k < channel_names.length; k++) {
					Table.set("Volume " + channel_names[k], row, volumes[k]);
				}				
			}
		}
	}	
	print("Number of chained clusters " + row);
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