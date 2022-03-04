#@File (label="folder", style="directory") top_folder
#@String(label="file extenstion",value="tif") file_extension
#@Boolean(label="Relative path",value=true) is_relative
/*
 * Parse sub folder and group images by subfolder names
 * 
 * Jerome Boulanger 2021
 */

condition_list = getFileList(top_folder);
n = nResults;
if (n != 0) {
	print("Appending new rows to existing table");
}
Array.print(condition_list);
for (c = 0; c < condition_list.length; c++) {
	print(condition_list[c]);
	condition_name = replace(condition_list[c],File.separator,'');
	image_list =  getFileList(top_folder + File.separator + condition_list[c]);
	Array.print(image_list);
	for (f = 0; f < image_list.length; f++) {
		if (endsWith(image_list[f], file_extension)) {			
			setResult("Condition",n,condition_name);
			if (is_relative) {
				setResult("Filename",n, condition_list[c] + image_list[f]);
			} else {
				setResult("Filename",n,top_folder + File.separator + condition_list[c] + image_list[f]);
			}
			n++;
		}
	}
}

//Table.rename('filelist.csv')
