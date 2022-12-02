#@String (label = "LUTs", description = "comma separated indices, or all",values="Red,Green,Blue") luts_csv
/*
 * Parse a sequence of text separated by comma
 * 
 * Jerome Boulanger
 */
 
default_luts  = newArray("Red","Green","Blue");
luts = parseCSVString(luts_csv, default_luts);
Array.print(luts);

function parseCSVString(csv, default) {	
	if (matches(csv, "auto")) {
		values = default;
	} else {
		str = split(csv,",");
		values = newArray(str.length);
		for (i = 0 ; i < str.length; i++) {
			values[i] = String.trim(str[i]);
		}
	}
	return values;
}