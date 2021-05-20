#@String(label="X Columns headers", value="Frame") xhead
#@String(label="Y Columns headers (csv or all)", value="ROI1-ch1,ROI1-ch2") yheadcsv
#@String(label="Y label", value="Intensity") ylabel
#@String(label="Palette", value="Ice") palette
#@String(label="Title", value="Graph") title
/*
 * Graph multi columns from the Result table
 * 
 * X Columns headers: the x serie
 * Y Columns headers: set to 'all' to plot all series, otherwise list the series (comma separated)
 * 
 * Jerome Boulanger 2021
 */

if (matches(yheadcsv,"all")) {
	yheadcsv = Table.headings;
	yheadtmp = split(yheadcsv,"\t");
	// remove empty and x header
	yhead = newArray();
	for (i=0;i<yheadtmp.length;i++) {
		if (!matches(yheadtmp[i],xhead) && !matches(yheadtmp[i]," ")) {
			yhead = Array.concat(yhead,yheadtmp[i]);
		}
	}
} else {
	yhead = split(yheadcsv,",");
}

Array.print(yhead);
showGraphs(xhead,yhead,palette,ylabel,title);

function showGraphs(xhead,yhead,palette,ylabel,title) {
	Plot.create(title,xhead,ylabel)
	legend="";
	colors = getLutHexCodes(palette,yhead.length);
	mini=255;
	maxi=0;
	T = Table.getColumn(xhead);
	for (c = 0; c < yhead.length; c++) {
		I = Table.getColumn(yhead[c]);
		Array.getStatistics(I, min, max, mean, stdDev);
		mini = minOf(min,mini);
		maxi = maxOf(max,maxi);
		legend += yhead[c]+"\t";
		Plot.setColor(colors[c]);
		Plot.setLineWidth(2);
		Plot.add("line", T, I,yhead[c]);
	}
	Plot.setLimits(0, T[T.length-1], mini, maxi);
	Plot.setLegend(legend,"options");
}

function getLutHexCodes(name,n) {
	setBatchMode(true);
	newImage("stamp", "8-bit White", 10, 10, 1);
   	run(name);
   	getLut(r, g, b);
   	close();
   	r = Array.resample(r,n);
   	g = Array.resample(g,n);
   	b = Array.resample(b,n);
	dst = newArray(r.length);
	for (i = 0; i < r.length; i++) {
		dst[i] = "#"+dec2hex(r[i]) + dec2hex(g[i]) + dec2hex(b[i]);
	}   
	setBatchMode(false);
	return dst;	
}

function dec2hex(n) {
	// convert a decimal number to its heximal representation
	codes = newArray("0","1","2","3","4","5","6","7","8","9","A","B","C","D","E","F");
	k1 = Math.floor(n/16);
	k2 = n - 16 * Math.floor(n/16);
	return codes[k1]+codes[k2]; 
}
