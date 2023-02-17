#@File("Input label image",style="open") inputFile

/*
 * Compute the rank distance to a point in the lanel image
 * Assume a clean label image 2D or 3D with only 1 channel and a ROI save in the TIFF file
 *
 * The reference label is defined by the current selection (ROI)
 * - Reports the shape measurement and rank information in a table
 * - Display an image of the labels by rank for visual inspection
 *
 * Jerome Boulanger for Nanami and Vivi 2023
 */

run("Close All");
setBatchMode("hide");
open(inputFile);
label_id = getImageID();

Roi.getPosition(foo, slice, foo);
setSlice(slice);
label = getValue("Mean");
print("Rank distance to label " + label);

ranks = rankToLabel(label);
rank_id = createRankImage(ranks);
rankfile = File.getNameWithoutExtension(inputFile) + "-ranks-to-label-"+label+".tif";
print("Saving " + rankfile);
saveAs("TIFF", File.getDirectory(inputFile) + File.separator + rankfile);
close();

tbl = createReport(label_id, ranks, label);
csvfile = File.getNameWithoutExtension(inputFile) + "-shape-and-ranks-to-label-"+label+".csv";
Table.rename(tbl, csvfile);
print("Saving " + csvfile);
Table.save(File.getDirectory(inputFile) + File.separator + csvfile);
selectWindow(csvfile);
run("Close");
close();
setBatchMode("exit and display");

function createReport(label_id, ranks, label) {
	/* Compute shape statistics on the label image and and the rank distance to the reference label
	 *  label_id (number) : id of the label image
	 *  ranks (array) : ranks for each label
	 */
	selectImage(label_id);
	is3d = (nSlices>1);
	if (is3d) {
		run("Analyze Regions 3D", "voxel_count volume surface_area mean_breadth sphericity euler_number bounding_box centroid equivalent_ellipsoid ellipsoid_elongations max._inscribed surface_area_method=[Crofton (13 dirs.)] euler_connectivity=6");
	} else {
		run("Analyze Regions", "area perimeter circularity euler_number bounding_box centroid equivalent_ellipse ellipse_elong. convexity max._feret oriented_box oriented_box_elong. geodesic tortuosity max._inscribed_disc average_thickness geodesic_elong.");
	}
	if (is3d) {
		line = 0;
		labels = Table.getColumn("Label");
		while (labels[line] != label && line < Table.size-1) {line++;}
		print("line:" + line + " label:" + labels[line]);
		X = Table.getColumn("Centroid.X");
		Y = Table.getColumn("Centroid.Y");
		Z = Table.getColumn("Centroid.Z");
	}
	for (i = 0; i < X.length; i++) {
		Table.set("Rank", i, ranks[i+1]);
		if (is3d) {
			dx = X[i] - X[line];
			dy = Y[i] - Y[line];
			dz = Z[i] - Z[line];
			Table.set("Euclidean Distance", i, sqrt(dx*dx+dy*dy+dz*dz));
		}
	}
	title = Table.title;
	return title;
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

	for (r = 1; r < rmax+1; r++) {
		list = "";
		for (i = 0; i < ranks.length; i++) {
			if (ranks[i] == r) {
				list += "" + (nlabels + i) + ",";
			}
		}
		if (list.length==0) break;
		list = substring(list, 0, list.length - 1);
		//print("" + r + ":" + list);
		run("Replace/Remove Label(s)", "label(s)="+list+" final="+r);
	}
	run("Replace/Remove Label(s)", "label(s)="+nlabels+" final=0");
	setMinAndMax(0, rmax);
	run("Fire");
	return getImageID();
}
