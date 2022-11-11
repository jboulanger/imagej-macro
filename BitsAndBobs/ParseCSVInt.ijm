#@String (label = "Angles", description = "comma separated indices, or all") angles_csv

default_angles  = newArray(1,2,3,4,5,6,7);
angles = parseCSVInt(angles_csv, default_angles);
Array.print(angles);

function parseCSVInt(csv, default) {
	print(csv);
	if (matches(csv, "all")) {
		values = default;
	} else {
		str = split(csv,",");
		values = newArray(str.length);
		for (i = 0 ; i < str.length; i++) {
			values[i] = parseInt(str[i]);
		}
	}
	return values;
}