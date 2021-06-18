#@ File (label="Input file", style="open") filename
#@ File (label="Output folder", style="directory") imageFolder
#@ String(label="Format",choices={"TIFF","PNG","JPEG"},style="list") format
#@ boolean (label="MIP",value=false) mip

/* Convert all series in a lif file to TIFs in a folder
 *  
 *  You can batch process LIFs using the Batch function in the script editor
 *  
 * Jerome Boulanger for Marta 2021
 */
 
run("Close All");
setBatchMode("hide");
run("Bio-Formats Macro Extensions");
print("File selected :" + filename);
Ext.setId(filename);
Ext.getSeriesCount(seriesCount);
print("File contains " + seriesCount + " series");
for (s = 1; s <= seriesCount; s++) {
	str="open=["+filename+"] color_mode=Composite rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_"+s+"";
	run("Bio-Formats Importer", str);
	if (mip) {
		print("Computing maximum intensity projection");
		id0=getImageID();
		run("Z Project...", "projection=[Max Intensity] all");	
		selectImage(id0);close();
	}
	oname = imageFolder + File.separator + replace(getTitle,".lif.*","_serie_"+IJ.pad(s,4)) + ".tif";
	print("Saving serie "+ s +"/" + seriesCount + " to "+ oname);
	setResult("",s-1,s)
	saveAs(format,oname);
	close();
}
Ext.close();
setBatchMode("exit and display");
