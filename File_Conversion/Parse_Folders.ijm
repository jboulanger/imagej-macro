#@File (label="folder", style="directory") top_folder
#@String(label="file extenstion",value="tif") file_extension
#@Boolean(lable="sub folder as group", value=false) group
#@Boolean(label="Relative path",value=true) is_relative
/*
 * Parse sub folder and group images by subfolder names
 *
 * Jerome Boulanger 2021
 */


if (group) {
	 parseGroupSubfolders(top_folder, file_extension, is_relative);
} else {
	parseFlatFolder(top_folder, file_extension, is_relative);
}

function parseFlatFolder(top_folder, file_extension, is_relative) {
	n = nResults;
	if (n != 0) {
		print("Appending new rows to existing table");
	}
	image_list =  getFileList(top_folder);
	Array.sort(image_list);
	for (f = 0; f < image_list.length; f++) {
		if (endsWith(image_list[f], file_extension)) {
			if (is_relative) {
				setResult("Filename",n, image_list[f]);
			} else {
				setResult("Filename",n,top_folder + File.separator + image_list[f]);
			}
			updateResults();
			n++;
		}
	}
	updateResults();
}

function parseGroupSubfolders(top_folder, file_extension, is_relative) {
	n = nResults;
	if (n != 0) {
		print("Appending new rows to existing table");
	}
	condition_list = getFileList(top_folder);
	for (c = 0; c < condition_list.length; c++) {
		condition_name = replace(condition_list[c],File.separator,'');
		image_list =  getFileList(top_folder + File.separator + condition_list[c]);
		for (f = 0; f < image_list.length; f++) {
			if (endsWith(image_list[f], file_extension)) {
				setResult("Condition",n,condition_name);
				if (is_relative) {
					setResult("Filename",n, condition_list[c] + image_list[f]);
				} else {
					setResult("Filename",n,top_folder + File.separator +
condition_list[c] + image_list[f]);
				}
				updateResults();
				n++;
			}
		}
	}
	updateResults();
}
