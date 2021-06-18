#@ File (label="Input file", style="open") filename
#@ File (label="Output folder", style="directory") imageFolder
#@ String(label="Format",choices={"TIFF","PNG","JPEG"},style="list") format
#@ boolean (label="MIP",value=false) mip
#@ boolean (label="Dummy run",value=true) dummy

/* Convert all series in a file to TIFs in a folder
 *  
 *  You can batch process files using the Batch function in the script editor
 *  
 *  Using the dummay mode enable to inspect image size as well and save the result in a cvs file
 *  
 * Jerome Boulanger for Marta 2021
 */
 
run("Close All");
setBatchMode("hide");
name = File.getNameWithoutExtension(filename);
ext = getNewFileExtension(format);


run("Bio-Formats Macro Extensions");
print("Input file:" + filename);
Ext.setId(filename);
Ext.getSeriesCount(seriesCount);
print("File contains " + seriesCount + " series");

for (s = 1; s <= seriesCount; s++) {
	
	str="open=["+filename+"] color_mode=Composite rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_"+s+"";	
	oname = imageFolder + File.separator + name + "_serie_" + IJ.pad(s,4) + ext;
	
	Ext.setSeries(s-1);
	Ext.getSeriesName(seriesName);
	Ext.getSizeX(sizeX);
	Ext.getSizeY(sizeY);
	Ext.getSizeZ(sizeZ);
	Ext.getSizeT(sizeT);
	if (dummy) {		
		print("bioformat args: " + str);	
		print("output file   : " + oname);		
		
	} else {
		
		print("Loading serie : " + s + "/" + seriesCount);			
		run("Bio-Formats Importer", str);
		
		if (mip) {
			print("Computing maximum intensity projection");
			id0 = getImageID();
			run("Z Project...", "projection=[Max Intensity] all");	
			selectImage(id0);close();
		}
		
		print("Saving serie "+ s +"/" + seriesCount + " to "+ oname);	
		saveAs(format,oname);
		close();		
	}
	setResult("Serie Index", s-1,s);
	setResult("Output  file", s-1, File.getName(oname));	
	setResult("Serie name", s-1, seriesName);
	setResult("Width", s-1, sizeX);
	setResult("Height", s-1, sizeY);
	setResult("Slices", s-1, sizeZ);
	setResult("Frames", s-1, sizeT);
	
}
Ext.close();
setBatchMode("exit and display");


function getNewFileExtension(format) {
	f = newArray("TIFF","PNG","JPEG");
	e = newArray(".tif",".png",".jpg");
	for (k = 0; k < f.length; k++) {
		if(matches(format, f[k])){
			k0 = k;
			return e[k];
		}
	}
	return ".tif";
}