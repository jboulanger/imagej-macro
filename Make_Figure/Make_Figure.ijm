//@String (label="Channel names", value="") channel_names_str
/*
 * Make a figure from a stack by splitting channels
 */
run("Close All");
if (nImages==0) {
 	loadTestImage();
}
channel_names = parse_channel_names(channel_names_str);
//setBatchMode("hide");
makeFigure(channel_names);
//setBatchMode("exit and display");

function makeFigure(channel_names) {
	id0 = getImageID();
	Stack.getDimensions(width, height, channels, slices, frames);
	getPixelSize(unit, pixelWidth, pixelHeight);
	print(slices);
	if (slices > 1) {
		run("Z Project...", "projection=[Max Intensity]");		
	} else {
		run("Duplicate...", "title=tmp duplicate");		
	}		
	id1 = getImageID();
	
	newImage("figure", "RGB black", (channels+1)*width, height, 1);
	id2 = getImageID();
	
	// create the composite overlay	
	selectImage(id1);
	Property.set("CompositeProjection", "Sum");
	Stack.setDisplayMode("composite");	
	run("Duplicate...", "title=tmp duplicate");		
	idtmp = getImageID();
	run("RGB Color");
	id3 = getImageID;
	if (idtmp!=id3) {selectImage(idtmp);close();selectImage(id3);}
	run("Copy");
	selectImage(id2);
	makeRectangle(width*(channels), 0, width, height);
	run("Paste");
	selectImage(id3);close();
	
	// copy each channel	
	selectImage(id1);
	Stack.setDisplayMode("grayscale");
	setColor("white");
	setFont("Arial", 25, "bold");
	
	for (c = 1; c <= channels; c++) {	
		selectImage(id1);
		Stack.setChannel(c);		
		run("Copy");
		selectImage(id2);
		makeRectangle(width*(c-1), 0, width, height);
		setPasteMode("Copy");
		run("Paste");
		if (c<channel_names.length-1) {
			Overlay.drawString(channel_names[c-1], width*(c-1)+5, 25);		
			Overlay.add();		
			Overlay.show();
		}
	}
	
	run("Select None");
	//selectImage(id1); close(); 
	
	Stack.setXUnit(unit);
	run("Properties...", "channels=1 slices=1 frames=1 pixel_width="+pixelWidth+" pixel_height="+pixelHeight+" voxel_depth=1");
	l = round(width * pixelWidth / 5);	
	run("Scale Bar...", "width="+l+" height=7 thickness=4 font=14 color=White background=None location=[Lower Right] horizontal bold overlay");
	return id2;
}

function parse_channel_names(channel_names_str) {
	channel_names = split(channel_names_str, ",");
	for (i = 0; i < channel_names.length; i++) { 
		channel_names[i] =  String.trim(channel_names[i]);
	}
	return channel_names;
}

function loadTestImage(){
	run("Clown");
	run("RGB Stack");
	run("Stack to Hyperstack...", "order=xyczt(default) channels=3 slices=1 frames=1 display=Color");
	Property.set("CompositeProjection", "Sum");	
	Stack.setDisplayMode("composite");
	Stack.setXUnit("cm");
	run("Properties...", "channels=3 slices=1 frames=1 pixel_width=0.1 pixel_height=0.1 voxel_depth=0.1");
}
