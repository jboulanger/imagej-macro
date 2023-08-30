/*
 * Create a list of objects
 * 
 * Set an object name
 * 
 */
 
var object_id = 0;

macro "set object id - C059T3e16ID" {	
	object_id = getNumber("Object ID", object_id + 1);
}

macro "add ROI [g] - C059T3e16+" {	
	roiManager("add");
	roiManager("rename", "obj" + IJ.pad(object_id, 4));
}