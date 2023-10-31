/*
 * This macro is called by Export cellSens VSI
 * It takes as argument a string x={xvalue} y={yvalue} giving
 * the position of the tile to process.
 * 
 * Results are appended in the result table.
 * 
 */
 
// parse the argument list
str = getArgument();
args = split(str);
for (i = 0; i < args.length; i++) {
	kv = split(args[i],"=");
	if (kv[0]=="x") {
		x = kv[1];
	} else if (kv[0]=="y") {
		y = kv[1];
	}
}

// process the tile (area over threshold in channel 1)
getDimensions(width, height, channels, slices, frames);
Stack.setChannel(1);
setAutoThreshold("Otsu");
run("Create Selection");
area = getValue("Area");
print(area);

row = nResults;
setResult("tile x", row, x);
setResult("tile y", row, y);
setResult("tile width", row, width);
setResult("tile height", row, height);
setResult("area above threshold", row, area);

// the macro is responsible for closing the tile
close();

