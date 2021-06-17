#@Integer(label="Columns", value=3) ncolumns
#@Integer(label="Rows", value=3) nrows
#@Integer(label="Horizontal overlap", value=8) overlap_x
#@Integer(label="Vertical overlap", value=8) overlap_y
#@String(label="Action", choices={"Split","Align","Merge"}) mode
/*
 * Split, Align, Merge tiles
 * 
 * 
 * Jerome Boulanger 2021
 */
run("Close All");
if (nImages==0) {
	print("testing mode");
	run("Close All");
	run("Blobs (25K)");
	//run("Boats");
	run("Grays");
	setBatchMode(true);
	split_tiles(ncolumns,nrows,overlap_x,overlap_y,10);
	setBatchMode(false);
	setBatchMode(true);
	align_tiles(ncolumns,nrows,overlap_x,overlap_y);
	setBatchMode(false);
	setBatchMode(true);
	merge_tiles(ncolumns,nrows,overlap_x,overlap_y);
	setBatchMode(false);
} else {
	setBatchMode(true);
	if (matches(mode,"Split")) {
		split_tiles(ncolumns,nrows,overlap_x,overlap_y,0);
	} else if (matches(mode,"Align")) {
		align_tiles(ncolumns,nrows,overlap_x,overlap_y);
	} else {
		merge_tiles(ncolumns,nrows,overlap_x,overlap_y);
	}
	setBatchMode(false);
}

function merge_tiles(ncolumns,nrows,overlap_x,overlap_y) {
	src = getImageID();
	Stack.getDimensions(tile_width, tile_height, channels, slices, frames);
	W = ncolumns*tile_width-(ncolumns-1)*overlap_x;
	H = nrows*tile_height-(nrows-1)*overlap_y;
	newImage("Tiled", "32-bit", W, H, 1);
	dst = getImageID();
	newImage("acc", "32-bit", W, H, 1);
	acc = getImageID();
	newImage("tmp1", "32-bit",  W, H, 1);
	tmp1 = getImageID();
	newImage("tmp2", "32-bit",  W, H,1);
	tmp2 = getImageID();
	b = 0;
	for (j = 0; j < nrows; j++) {
		for (i = 0; i < ncolumns; i++) {
			x =  i * (tile_width - overlap_x);
			y =  j * (tile_height - overlap_y);
			selectImage(src);
			Stack.setFrame(i+ncolumns*j+1);
			makeRectangle(b,b,tile_width-2*b,tile_height-2*b);
			run("Copy");
					
			// create a mask for the tile
			selectImage(tmp2);
			run("Select None");
			setColor(0); fill();
			makeRectangle(x+b,y+b,tile_width-2*b,tile_height-2*b);
			run("Paste");
			run("Select None");
			run("Macro...", "code=v=(v>1) slice");
			// add the tmp image to the accumulator acc
			imageCalculator("Add 32-bit", acc, tmp2);

			// paste the tile in tmp1
			selectImage(tmp1);
			run("Select None");
			setColor(0); fill();
			makeRectangle(x+b,y+b,tile_width-2*b,tile_height-2*b);
			run("Paste");
			run("Select None");
			imageCalculator("Multiply 32-bit", tmp1, tmp2);
			// add the tmp image to the accumulator dst
			imageCalculator("Add", dst, tmp1);
		}
	}
	imageCalculator("Divide", dst, acc);
	selectImage(acc); close();
	selectImage(tmp1); close();
	selectImage(tmp2); close();
	selectImage(dst);
	run("Select None");
	return dst;
}

function align_tiles(ncolumns,nrows,overlap_x,overlap_y) {
	src = getImageID();
	Stack.getDimensions(tile_width, tile_height, channels, slices, frames);
	ref = newArray(nrows*ncolumns);
	idx = newArray(nrows*ncolumns);
	sideref = newArray(nrows*ncolumns);
	sideidx = newArray(nrows*ncolumns);
	createPath(ref,idx,sideref,sideidx,ncolumns,nrows,0);
	processPath(ref,idx,sideref,sideidx,overlap_x,overlap_y);
	createPath(ref,idx,sideref,sideidx,ncolumns,nrows,2);
	processPath(ref,idx,sideref,sideidx,overlap_x,overlap_y);
	//createPath(ref,idx,sideref,sideidx,ncolumns,nrows,1);
	//processPath(ref,idx,sideref,sideidx,overlap_x,overlap_y);
	run("Select None");
}

function processPath(ref,idx,sideref,sideidx,overlap_x,overlap_y) {
	// register images following the path defined by the arguments
	for (k = 1; k < ref.length; k++) {
		shift = cropAndEstimateShift(ref[k],idx[k],sideref[k],sideidx[k],overlap_x,overlap_y);
		Stack.setFrame(idx[k]);
		run("Select None");
		run("Translate...", "x="+-shift[0]+" y="+-shift[1]+" interpolation=Bicubic");
	}
}

function createPath(ref,idx,sideref,sideidx,ncolumns,nrows,type) {
	// create a path across tiles with matching sides
	k = 0;
	if (type==0) {
		for (j = 0; j < nrows; j++) {
			for (i = 0; i < ncolumns; i++) {
				if (k!=0) {
					idx[k] = i+ncolumns*j+1;
					if (j==0) {
						ref[k] = (i-1)+ncolumns*j+1;
						sideref[k] = 2;
						sideidx[k] = 4;
					} else {
						ref[k] = i+ncolumns*(j-1)+1;
						sideref[k] = 3;
						sideidx[k] = 1;
					}
				}
				k++;
			}
		}
	} else if(type==1) {
		for (i = 0; i < ncolumns; i++) {
			for (j = 0; j < nrows; j++) {	
				if (k!=0) {
					idx[k] = i+ncolumns*j+1;
					if (i==0) {
						ref[k] = i+ncolumns*(j-1)+1;
						sideref[k] = 3;
						sideidx[k] = 1;
					} else {
						ref[k] = i-1+ncolumns*j+1;
						sideref[k] = 2;
						sideidx[k] = 4;
					}
				}
				k++;
			}
		}
	} else if (type==2) {
		for (j = nrows-1; j > 0; j--) {
			for (i = ncolumns-1; i > 0; i--) {
				if (k!=nrows*ncolumns-1) {
					idx[k] = i+ncolumns*j+1;
					if (j==nrows-1) {
						ref[k] = (i+1)+ncolumns*j+1;
						sideref[k] = 4;
						sideidx[k] = 2;
					} else {
						ref[k] = i+ncolumns*(j+1)+1;
						sideref[k] = 1;
						sideidx[k] = 3;
					}
				}
				k++;
			}
		}
	}
}

function estimateShift(id1,id2,scale) {
	sx = 0;
	sy = 0;
	n = 1;
	for (iter = 0; iter < 20 && n > 0.01; iter++) {
		selectImage(id2); run("Duplicate...", " "); id3 = getImageID();
		run("Translate...", "x="+-sx+" y="+-sy+" interpolation=Bilinear");
		shift = estimateShiftCore(id1,id3,scale);
		selectImage(id3); close();
		sx -= shift[0];
		sy -= shift[1];
		n = sqrt(shift[0]*shift[0]+shift[1]*shift[1]);
	}
	return newArray(sx,sy);
}

function estimateShiftCore(id1,id2,scale) {
	selectImage(id1);
	run("Duplicate...", " ");
	run("32-bit");
	run("Gaussian Blur...", "sigma="+scale);
	ix = getImageID();rename("ix");
	selectImage(id2);
	run("Duplicate...", " ");
	run("32-bit");
	run("Gaussian Blur...", "sigma="+scale);
	iy = getImageID();rename("iy");

	imageCalculator("Subtract create 32-bit", id2, id1);
	it = getImageID();rename("it");
	
	selectImage(ix);
	run("Convolve...", "text1=[0 0 0\n-1 0 1\n0 0 0\n]");
	
	selectImage(iy);
	run("Convolve...", "text1=[0 -1 0\n0 0 0\n0 1 0\n]");
	// compute the second order images ixx iyy 
	imageCalculator("Multiply create 32-bit", ix, it); ixt = getImageID();rename("ixt");
	imageCalculator("Multiply create 32-bit", iy, it); iyt = getImageID();rename("iyt");
	imageCalculator("Multiply create 32-bit", ix, iy); ixy = getImageID();rename("ixy");
	selectImage(ix); run("Square");
	selectImage(iy); run("Square");
	// create a mask
	imageCalculator("Multiply create 32-bit", id1, id2); mask = getImageID();rename("mask");
	run("Macro...", "code=v=(v>1) slice");
	ids = newArray(ix,iy,it,ixy,ixt,iyt);
	for (i = 0; i < ids.length; i++) {imageCalculator("Multiply", ids[i], mask);}
	// measure the mean intensity in each image
	getDimensions(w, h, channels, slices, frames);
	band = 2;
	selectImage(ix);  makeRectangle(band, band, w-2*band, h-2*band); a = getValue("Mean");
	selectImage(ixy); makeRectangle(band, band, w-2*band, h-2*band); b = getValue("Mean");
	selectImage(iy);  makeRectangle(band, band, w-2*band, h-2*band); c = getValue("Mean");
	selectImage(ixt); makeRectangle(band, band, w-2*band, h-2*band); d = getValue("Mean");
	selectImage(iyt); makeRectangle(band, band, w-2*band, h-2*band); e = getValue("Mean");
	// solve [a,b;b,c] [vx;vy] = [d;e] 
	vx = 2*(c*d-b*e) / (a*c-b*b);
	vy = 2*(a*e-b*d) / (a*c-b*b);
	// clear the list of images
	for (i = 0; i < ids.length; i++) {selectImage(ids[i]);close();}
	selectImage(mask);close();
	return newArray(vx,vy);
}

function cropAndEstimateShift(n1,n2,side1,side2,overlap_x,overlap_y) {
	id0=getImageID();
	run("Select None");
	id1 = duplicate_side("A",n1,side1,overlap_x,overlap_y);
	selectImage(id0);
	id2 = duplicate_side("B",n2,side2,overlap_x,overlap_y);

	shift = estimateShift(id1,id2,0);
	selectImage(id1); close();
	selectImage(id2); close();
	return shift;
}


function duplicate_side(name,n,side,overlap_x,overlap_y) {
	Stack.getDimensions(width, height, channels, slices, frames);
	Stack.setFrame(n);
	if (side==1) {
		makeRectangle(0, 0, width, overlap_y);
	} else if (side==2) {
		makeRectangle(width-overlap_x-1, 0, overlap_y, height);
	} else if (side==3) {
		makeRectangle(0, height-overlap_y-1, width, overlap_y);
	} else if (side==4) {
		makeRectangle(0, 0, overlap_x, height);
	} else {
		run("Select None");
	}
	run("Duplicate...", "title="+name);
	run("Unsharp Mask...", "radius=2 mask=0.9");
	return getImageID();
}

function split_tiles(ncolumns,nrows,overlap_x,overlap_y,noise) {
	src = getImageID();
	Overlay.remove();
	Stack.getDimensions(width, height, channels, slices, frames);
	tile_width = Math.floor((width+(ncolumns-1)*overlap_x) / ncolumns);
	tile_height = Math.floor((height+(nrows-1)*overlap_y) / nrows);
	newImage("Tiles", "8-bit", tile_width, tile_height, channels, slices, ncolumns*nrows);
	dst = getImageID();
	for (j = 0; j < nrows; j++) {
		for (i = 0; i < ncolumns; i++) {
			x = i * (tile_width - overlap_x) + noise*random;
			y = j * (tile_height - overlap_y) + noise*random;
			for (channel = 1; channel <= channels; channel++) {
				for (slice = 1; slice <= slices; slice++) {
					selectImage(src);
					Stack.setPosition(channel, slice, 1);
					makeRectangle(x,y,tile_width,tile_height);
					Stack.setPosition(channel, slice, 1);
					run("Copy");
					selectImage(dst);
					Stack.setPosition(channel,slice,i+ncolumns*j+1);
					run("Paste");
				}
			}
		}
	}
	run("Select None");
}

