#@Integer (label="cell channel",value=1) chan1
#@Integer (label="nuclei channel",value=2) chan2
#@Float (label="diameter", value=30) diameter
#@String(label="model", choices={"cyto2","LC2"}) model
#@Boolean(label="use GPU", value=true) use_gpu
#@File (label="temporary folder", style="directory") tmppath
#@String(label="emvironment", value="cellpose") environment

/*
 * Run cellpose on the current image
 */


cellpose2d(tmppath,chan1,chan2,model,diameter,use_gpu,environment);

function  cellpose2d(tmpdir,ch1,ch2,model,diameter,use_gpu,environment) {
	ofile = tmpdir + File.separator + "tmp.tif";
	mfile = tmpdir + File.separator + "tmp_cp_outlines.txt";
	
	if (ch2==0) {
		run("Duplicate...", "duplicate channels=" + ch1);
	} else {
		run("Make Subset...", "channels="+ch1+","+ch2);
	}
		
	print("Saving tiff file for segmentation\n" + ofile);
	saveAs("TIFF", ofile);
	close();
	cmd = "conda run -n "+environment+" python -m cellpose --chan " + ch1 + " --chan2 "+ch2;
	cmd += " --image_path " + ofile;
	cmd += " --pretrained_model "+model+" --diameter " + diameter;
	cmd += " --save_txt";
	cmd += " --use_gpu";
	print("running cellpose command " + cmd);
	ret = exec(cmd);
	print("Anszwer:"+ret);
	print("loading outlines");
	loadCellposeOutlines(mfile);
	File.delete(ofile);
	File.delete(mfile);
}


function loadCellposeOutlines(path) {
	str = File.openAsRawString(path);	
	listcoords = split(str, "\n");			
	for (i = 0; i < listcoords.length; i++) {
		xy = split(listcoords[i],',');
		n = xy.length/2;
		x = newArray(n);
		y = newArray(n);
		for (j = 0; j < n; j++) {
			x[j] = xy[2*j];
			y[j] = xy[2*j+1];
		}
		makeSelection("freehand", x, y);
		roiManager("add");		
	}
	roiManager("Show All");
}
