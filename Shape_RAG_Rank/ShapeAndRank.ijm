#@File("Input label image",style="open") inputFile

/* 
 * Compute the rank distance to a point in the lanel image
 * Assume a clean label image
 * 
 * The reference label is defined by the current selection (ROI)
 * - Reports the shape measurement and rank information in a table
 * - Display an image of the labels by rank for visual inspection
 * 
 * Jerome Boulanger for Nanami and Vivi 2023
 */
 
open(inputFile);

label_id = getImageID();
label = getValue("Mean");

print("Rank distance to label " + label);
ranks = rankToLabel(label);
rank_id = createRankImage(ranks);
rankfile = File.getNameWithoutExtension(inputFile) "-ranks.tif";
saveAs(rankfile);

createReport(label_id, ranks);
csvfile = File.getNameWithoutExtension(inputFile) "-shape-and-ranks.csv";
saveAs(csvfile);

function createReport(label_id, ranks, label) {
	/* Compute shape statistics on the label image and and the rank distance to the reference label
	 *  label_id (number) : id of the label image
	 *  ranks (array) : ranks for each label
	 */
	selectImage(label_id);
	if (nSlices > 1) {
		run("Analyze Regions 3D", "voxel_count volume surface_area mean_breadth sphericity euler_number bounding_box centroid equivalent_ellipsoid ellipsoid_elongations max._inscribed surface_area_method=[Crofton (13 dirs.)] euler_connectivity=6");
	} else {
		run("Analyze Regions", "area perimeter circularity euler_number bounding_box centroid equivalent_ellipse ellipse_elong. convexity max._feret oriented_box oriented_box_elong. geodesic tortuosity max._inscribed_disc average_thickness geodesic_elong.");
	}
	for (i = 0; i < Table.size; i++) {
		Table.set("Rank", i, ranks[i+1]);
	}
	Table.rename(name + "rank-distance-to-label-" + label +".csv");	
}

function rankToLabel(index) {
	/* Compute the rank distance to a label
	 *  index (number) : label
	 * Returns (array) : rank for each label
	 */
	print("Computing rank distance");
	run("Select None");
	if (nSlices > 1) {		
		Stack.getStatistics(voxelCount, mean, min, max, stdDev);
	} else {		
		getStatistics(area, mean, min, max, std, histogram);
	}
	nlabels = max + 1;
	
	run("Region Adjacency Graph", " image=["+getTitle()+"]");

	a = Table.getColumn("Label 1");
	b = Table.getColumn("Label 2");
	n = a.length;
		
	// Close the table with the RAG graph
	run("Close");
	
	ranks = newArray(nlabels);

	ranks[index] = 1;

	for (r = 2; r < 100; r++) {
		tmp = newArray(nlabels);  // list of neightbords
		k = 0;                    // number of neightbors
		for (i = 0; i < n; i++) { // for each edge of the RAG 
			ra = ranks[a[i]];     // get the ranks of each vertex
			rb = ranks[b[i]];
			if ( ra > 0 && rb == 0) { // add bi if on the border
				tmp[k] = b[i];
				k++;
			} else if (rb > 0  && ra == 0) { // add ai if on the border
				tmp[k] = a[i];
				k++;
			}
		}
		if (k==0) {break;}
		for (i = 0; i < k; i++) {
			ranks[tmp[i]] = r;
		}
	}
	return ranks;
}

function createRankImage(ranks) {
	/* Create an image with labeled as ranks */
	print("Create a rank distance image");
	nlabels = ranks.length;
	Array.getStatistics(ranks, min, rmax, mean, stdDev);
	run("Duplicate...", "title=rank duplicate");
	run("Add...", "value=" + nlabels + " stack");
	Array.print(ranks);
	for (r = 1; r < rmax+1; r++) {
		list = "";
		for (i = 0; i < ranks.length; i++) {
			if (ranks[i] == r) {
				list += "" + (nlabels + i) + ",";
			}
		}		
		if (list.length==0) break;		
		list = substring(list, 0, list.length - 1);
		print("" + r + ":" + list);
		run("Replace/Remove Label(s)", "label(s)="+list+" final="+r);
	}	
	run("Replace/Remove Label(s)", "label(s)="+nlabels+" final=0");			
	setMinAndMax(0, rmax);	
	run("Fire");
	return getImageID();
}


