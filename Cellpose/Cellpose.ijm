//@Integer (label="cell channel",value=1) chan1
//@Integer (label="nuclei channel",value=2) chan2
//@Float (label="diameter", value=30) diameter
//@String(label="model", choices={"cyto2","cyto3","LC2"}) model
//@Boolean(label="use GPU", value=true) use_gpu

/*
 * Run cellpose on the current image
 * 
 * Jerome Boulanger 2023
 */

/* Configure here cellpose */
conda = "micromamba";
// conda = "C:\\Users\\User\\miniconda3\Script\conda.exe"
environment = "imaging";
tmpdir = "/tmp";
use_gpu = true;
/* end of configuration */

if (nImages==0) {	
	print("No image opened, loading Fluorescent Cells example with automatic settings");
	run("Fluorescent Cells");	
	chan1 = 1;
	chan2 = 3;
	diameter = 100;	
}

cellpose2d(chan1,chan2,model,diameter);

function  cellpose2d(ch1,ch2,model,diameter) {
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
	
	if (ch2==0) {
		cmd = conda + " run -n "+environment+" python -m cellpose --chan 1 --chan2 0";
	} else {
		cmd = conda + " run -n "+environment+" python -m cellpose --chan 1 --chan2 2";
	}	
	cmd += " --image_path " + ofile;
	cmd += " --pretrained_model "+model+" --diameter " + diameter;
	cmd += " --save_txt";
	cmd += " --use_gpu";
	print("running cellpose command " + cmd);
	ret = exec(cmd);
	print("Anszwer:"+ret);
	print("loading outlines from file "+ mfile);
	loadCellposeOutlines(mfile);
	File.delete(ofile);
	File.delete(mfile);
}


function loadCellposeOutlines(path) {
	if (File.exists(path)) {
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
			Roi.setStrokeColor("yellow");
			roiManager("add");		
		}
		roiManager("Show All");
	} else {
		print("No result file found");
	}
}
