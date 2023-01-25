#@File (label="Input file", style="open") filename
#@Float (label="Zoom factor", value=0.25) zoom
/*
 * Tile a multiposition file using the metadata.
 * 
 * Jerome Boulanger 2022
 */
 
run("Bio-Formats Macro Extensions");
print("Input file:" + filename);
Ext.setId(filename);
Ext.getSeriesCount(seriesCount);

print("File contains " + seriesCount + " series");
X = newArray(seriesCount);
Y = newArray(seriesCount);
Z = newArray(seriesCount);
W = newArray(seriesCount);
H = newArray(seriesCount);
D = newArray(seriesCount);
x0 = 1e9;
y0 = 1e9
z0 = 1e9;
for (s = 0; s < seriesCount; s++) {		
	Ext.setSeries(s);
	Ext.getSizeX(nx);
	Ext.getPixelsPhysicalSizeX(dx);
	Ext.getSizeY(ny);
	Ext.getPixelsPhysicalSizeY(dy);
	Ext.getSizeZ(nz);
	Ext.getPixelsPhysicalSizeZ(dz);
	Ext.getSizeC(channels);
	Ext.getSeriesMetadataValue("Image #0|Tile #"+s+"|PosX",x);
	Ext.getSeriesMetadataValue("Image #0|Tile #"+s+"|PosY",y);
	Ext.getSeriesMetadataValue("Image #0|Tile #"+s+"|PosZ",z);	
	X[s] = 1000000 * x / dx;
	Y[s] = 1000000 * y / dy;
	Z[s] = 1000000 * z / dz;
	W[s] = nx;
	H[s] = ny;
	D[s] = nz;
	x0  = minOf(X[s], x0);
	y0  = minOf(Y[s], y0);
	z0  = minOf(Z[s], z0);
}
width = 0;
height = 0;
depth = 0;
for (s = 0; s < seriesCount; s++) {		
	X[s] = X[s] - x0;
	Y[s] = Y[s] - y0;
	Z[s] = Z[s] - z0;
	width = maxOf(X[s]+W[s], width);	
	height = maxOf(Y[s]+H[s], height);	
	depth = maxOf(Z[s]+D[s], depth);	
}

width = Math.ceil(width);
height = Math.ceil(height);
depth = Math.ceil(depth);
Array.show("positions",X,Y,Z,W,H,D);
print("Image size " + width + " x " + height + " x " + depth + " x " + channels);

// mosaic
setBatchMode("hide");
print("Image size " + width + " x " + height + " x " + depth + " x " + channels);
if (depth > 1/zoom) {	
	newImage("tiled", "16-bit", round(zoom*width) , round(zoom*height), channels, round(zoom*depth), 1);
} else {
	newImage("tiled", "16-bit", zoom*width , zoom*height, channels, depth, 1);
}
setPasteMode("Blend");
id0 = getImageID();
for (s = 0; s < seriesCount; s++) {
	Ext.setSeries(s);
	w = round(zoom * W[s]);
	h = round(zoom * W[s]);
	d = round(zoom * D[s]);
	dz = round(1./zoom);
	for (z = 0; z < D[s]; z = z + dz) {
		for (c = 0; c < channels; c++) {	
			showProgress(c+channels*(z+zoom*depth*s), channels*zoom*depth*seriesCount);
			print("serie " + s + " plane " + z + " channel " + c);
			Ext.openImage("tmp", c+z*channels);
			info = getMetadata("Info");
			run("Size...", "width="+w+" height="+h+" depth=1 constrain average interpolation=None");					
			run("Rotate 90 Degrees Right");			
			run("Copy");
			wcrop = getWidth();
			hcrop = getHeight();
			close();
			selectImage(id0);
			Stack.setPosition(c+1, round(zoom * z), 1);
			makeRectangle(round(zoom*X[s]), round(zoom*Y[s]), wcrop, hcrop);
			run("Paste");
		}
	}
	Overlay.drawRect(zoom*X[s], zoom*Y[s], w, d);
	Overlay.setPosition(0, 0, 0);
	Overlay.add();
}
Overlay.show();
setBatchMode("exit and display");
setMetadata("Info", info);
