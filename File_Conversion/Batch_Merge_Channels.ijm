#@File(label="Filename") filename
#@String(Label="Channel template",description="comma separated list of pattern to replace for guessing the channel file",value="C1,C2,C3,C4") channel_template_str
#@String(Label="Output file tag",description="tab inserted before the file extension when saving the result",value="-merged") tag
/*
 * Merge channels from a set of files
 * 
 * Replace template string of the 1st listed channel by the subsequent one
 * and merge all opened images into a single one. Save the result.
 *
 * Parameters:
 * Filename : Select the file corresponding to the first channel eg: C1-mypicture.tif (or if you selected batch for running the macro all the files corresponding to the 1st channel)
 * Channel template: list with a comma without space, the part of the file name which match the first one: eg C1,C2,C3
 * Output file tag: a string to add to the filename when saving the result (-merged.tif)
 * 
 * Jerome for Giulia 2022
 */
 
channel_template = split(channel_template_str,",");
n = channel_template.length;
str = "";

setBatchMode(true);
for (k = 0; k < n; k++) {	
	f = replace(filename,channel_template[0],channel_template[k]);	
	open(f);	
	t = getTitle();
	str += "c"+(k+1) +"=["+t+"]";
}

run("Merge Channels...", str + " create");

opath =  File.getDirectory(filename) + File.separator + File.getNameWithoutExtension(filename) + tag + ".tif";
print("Saving file \n" + opath);
saveAs("TIFF", opath);
close();
setBatchMode(false);