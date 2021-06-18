#@ File (label="Input file", style="open") filename
#@ File (label="Output folder", style="directory") imageFolder

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
for (s = 1; s <= seriesCount; s++) {
	str="open=["+filename+"] color_mode=Composite rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_"+s+"";
	run("Bio-Formats Importer", str);
	oname = imageFolder + File.separator + replace(getTitle,".lif","_serie_"+s) + ".tif";
	print("Saving serie "+ s +"/" + seriesCount + " to "+ oname);
	saveAs("TIFF",oname);
	close();
}
Ext.close();
setBatchMode("exit and display");
