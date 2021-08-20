#@File (label="folder", style="directory") top_folder
#@String(label="file extenstion",value=".lsm") file_extension
/*
 * Parse sub folder and group images by subfolder names
 * 
 * 
 */

condition_list = getFileList(top_folder);
n = nResults;
if (n != 0) {
	print("Appending new rows to existing table");
}
for (c = 0; c < condition_list.length; c++) {
	condition_name = replace(condition_list[c],File.separator,'');
	image_list =  getFileList(top_folder + File.separator + condition_list[c]);
	for (f = 0; f < image_list.length; f++) {
		if (endsWith(image_list[f], file_extension)) {
			setResult("condition",n,condition_name);
			setResult("filename",n,top_folder + File.separator + condition_list[c] + image_list[f]);
			n++;
		}
	}
}


