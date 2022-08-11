/*
 * Measure sum of the intensity along time in various ROI and compute the ratio 
 * of each versus the sum of the others.
 * 
 * Jerome Boulanger 2022 for Yangci
 */
 
start_time = getTime();
tbl = "Intensity fraction.csv";
//setBatchMode("hide");
measureAllROIsIntensitiesOverTime(tbl);
plotFracVsTime(tbl);
//setBatchMode("exit and display");
end_time = getTime();
print("Elapsed time " + (end_time - start_time)/1000 + "s");
 
  
function measureAllROIsIntensitiesOverTime(tbl) {
	n = roiManager("count");
	Stack.getDimensions(width, height, channels, slices, frames);
	tunits = getTimeUnits();
	Table.create(tbl);		
	for (frame = 1; frame <= frames; frame++) {
		Table.set(tunits[2], frame-1, (frame-1) * tunits[0]);
		for (channel = 1; channel <= channels; channel++) {
			Stack.setPosition(channel, 1, frame);
			sum = 0;
			vals = newArray(n);
			for (index = 0; index < n; index++) {
				roiManager("select", index);				
				Stack.setPosition(channel, 1, frame);
				vals[index] = getValue("RawIntDen");
				sum += vals[index];				
			}
			for (index = 0; index < n; index++) {
				cname = "Frac.R"+(index+1)+"C"+channel;
				Table.set(cname, frame-1, vals[index]/sum);
			}			
		}
	}
}

function getTimeUnits() {
	// return an array with [ FrameInterval,Unit,"Time ["+Unit"]" ]
	dt = Stack.getFrameInterval();
	if (dt==0) {
		dt = 1;
		Time = "Frame";
	} else {
		Stack.getUnits(X, Y, Z, Time, Value);
	}
	return newArray(dt, Time, "Time ["+Time+"]");
}

function plotFracVsTime(tbl) {	
	// plot all intensity fraction versus time from the table
	n = roiManager("count");
	Stack.getDimensions(width, height, channels, slices, frames);	
	tunits = getTimeUnits();
	time = Table.getColumn(tunits[2]);
	colors = getAllChannelsLutTriplet();
	Array.print(colors);
	Plot.create("Fractions", tunits[2], "intensity fractions");
	Plot.setLineWidth(2);
	legend="";
	for (channel = 1; channel <= channels; channel++) {		
		for (index = 0; index < n; index++) {	
			cname = "Frac.R"+(index+1)+"C"+channel;
			values = Table.getColumn(cname);	
			col = createColorHexCode(colors[channel-1], index, n, 0.45, 0.75);
			Plot.setColor(col);			
			Plot.add("line",time,values,cname);
			legend += cname+"\t";			
		}
	}
	Plot.setLegend(legend, "bottom-right");
	Plot.setLimits(0, time[time.length-1], 0, 1);
	Plot.update();	
}

function dec2hex(n) {
	// convert a decimal number to its heximal representation
	codes = newArray("0","1","2","3","4","5","6","7","8","9","A","B","C","D","E","F");
	k1 = Math.floor(n/16);
	k2 = n - 16 * Math.floor(n/16);
	return codes[k1]+codes[k2]; 
}

function getAllChannelsLutTriplet() {
	// compute array with colors as hex code for all channels 
	Stack.getDimensions(width, height, channels, slices, frames);
	colors = newArray(channels);
	for (channel = 1; channel <= channels; channel++) {
		Stack.setChannel(channel);
		getLut(reds, greens, blues);
		n = reds.length-1;		
		colors[channel-1] = ""+reds[n]+","+greens[n]+","+blues[n];
	}
	return colors;
}

function createColorHexCode(color_triplet, index, n, low, high) {
	// create a color gradient based on index
	c = split(color_triplet,",");
	Array.print(c);
	q = low + (high-low) * index / (n-1);	
	code = "#";
	for (k = 0; k < c.length; k++) {
		x = round(q*c[k]);
		print(x);
		code = code + dec2hex(round(q*c[k]));
	}
	return code;	
}

/*
function getChannelColorAsHex() {
	// get the lut of a channel and return the brightest color as hex color code
	getLut(reds, greens, blues);
	n = round(0.75*reds.length);	
	return"#"+dec2hex(reds[n])+dec2hex(greens[n])+dec2hex(blues[n]);	
}

function getAllChannelColorsAsHex() {
	// compute array with colors as hex code for all channels 
	Stack.getDimensions(width, height, channels, slices, frames);
	colors = newArray(channels);
	for (channel = 1; channel <= channels; channel++) {
		Stack.setChannel(channel);
		str = getChannelColorAsHex();
		colors[channel-1] = str;
	}
	return colors;
}
*/
