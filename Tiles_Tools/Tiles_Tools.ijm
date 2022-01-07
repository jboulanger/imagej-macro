#@Integer(label="Columns", value=3) ncolumns
#@Integer(label="Rows", value=3) nrows
#@Integer(label="Horizontal overlap", value=8) overlap_x
#@Integer(label="Vertical overlap", value=8) overlap_y
#@Integer(label="Position noise", value=0) noise
#@String(label="Action", choices={"Split","Distribute","Update","Align","Apply","Correct","Merge"}) mode

/*
 * Split, Align & Merge tiles
 * 
 * Jerome Boulanger 2021
 */
 
if (nImages==0) {
	print("testing mode");
	run("Close All");
	run("Blobs (25K)");
	//run("Boats");
	run("Grays");
	setBatchMode(true);
	print("Splitting the image into overlaping tiles");
	tiles = split_tiles(ncolumns,nrows,overlap_x,overlap_y,2);
	setBatchMode(false);
		
	print("Distributing tiles");
	view = spread_tiles(ncolumns,nrows,overlap_x,overlap_y);
	close();	
	setBatchMode(true);
	print("Aligning the images tiles");
	selectImage(tiles);
	align_tiles(ncolumns,nrows,overlap_x,overlap_y);
	setBatchMode(false);
	
	setBatchMode(true);
	print("Merging the tiles in a mosaic");
	merge_tiles(ncolumns,nrows,overlap_x,overlap_y);
	setBatchMode(false);
} else {
	setBatchMode("hide");
	if (matches(mode,"Split")) {
		split_tiles(ncolumns,nrows,overlap_x,overlap_y,noise);
	} else if (matches(mode,"Align")) {
		nalign_tiles(ncolumns,nrows,overlap_x,overlap_y);
	} else if (matches(mode,"Apply")) {
		apply_transform(ncolumns,nrows,overlap_x,overlap_y);
	} else if (matches(mode,"Merge")){
		merge_tiles(ncolumns,nrows,overlap_x,overlap_y);
	} else if (matches(mode,"Distribute")) {
		spread_tiles(ncolumns,nrows,overlap_x,overlap_y);
	} else if (matches(mode,"Update")) {
		readout(ncolumns,nrows,overlap_x,overlap_y);
	} else if (matches(mode,"Correct")) {
		correct_illumination();
	}	
	setBatchMode("exit and display");
}

function spread_tiles(ncolumns,nrows,overlap_x,overlap_y) {	
	// spread the tiles from a stack
	setBatchMode("exit and display");
	if (isOpen("ROI Manager")) {selectWindow("ROI Manager");run("Close");}
	id0 = getImageID;
	name = getTitle();
	Stack.getDimensions(tile_width, tile_height, channels, slices, frames);	
	W = round(1.1*ncolumns*tile_width-(ncolumns-1)*overlap_x);
	H = round(1.1*nrows*tile_height-(nrows-1)*overlap_y);
	newImage("View", "8-bit composite-mode", W, H, channels, slices, 1);
	id1 = getImageID();	
	for (j = 0; j < nrows; j++) {
		for (i = 0; i < ncolumns; i++) {
			k = i + nrows * j;
			x = round(i * (tile_width - overlap_x) + 0.05 * ncolumns * tile_width);
			y = round(j * (tile_height - overlap_y) + 0.05 * nrows * tile_height);
			selectImage(id0);
			Stack.setPosition(1, 1, k + 1);
			selectImage(id1);			
			run("Add Image...", "image=["+name+"] x="+x+" y="+y+" opacity=50");					
		}
	}	
	run("To ROI Manager");
	readout(ncolumns,nrows,overlap_x,overlap_y);
	return id1;
}

function getTileCoords() {
	// get the tiles coordinates
	return newArray(getValue("BX"), getValue("BY"), getValue("Width"), getValue("Height"));
}

function getTileCoordsFromTable(row) {
	// get the tiles coordinates
	if (row < nResults) {
		return newArray(getResult("x", row),getResult("y", row),getResult("width", row),getResult("height", row));
	} else {
		print("getTileCoordsFromTable row "+ row + " out of range " + nResults);
		return newArray(0,0,0,0);
	}
}

function areaTilesOverlap(A,B) {
	// return the area of the overlap (intersection) between the two tiles
	return maxOf(0, minOf(A[0]+A[2], B[0]+B[2])-maxOf(A[0], B[0])) * maxOf(0, minOf(A[1]+A[3], B[1]+B[3])-maxOf(A[1], B[1]));	
}

function readout(ncolumns,nrows,overlap_x,overlap_y) {
	// read the position of the tiles and populate to the result table
	n = roiManager("count");
	for (j = 0; j < nrows; j++) {
		for (i = 0; i < ncolumns; i++) {
			k = i + nrows * j;
			roiManager("select", k);						
			x = getValue("BX");
			y = getValue("BY");
			w = getValue("Width");
			h = getValue("Height");
			setResult("ROI", k, Roi.getName);
			setResult("Position", k, k + 1);
			setResult("x", k, x);
			setResult("y", k, y);
			setResult("width", k, w);
			setResult("height", k, h);
			setResult("dx", k, 0);
			setResult("dy", k, 0);
		}
	}	
	roiManager("Show All");	
}

function merge_tiles(ncolumns,nrows,overlap_x,overlap_y) {
	// merge tiles in a mosaic
	src = getImageID();
	Stack.getDimensions(tile_width, tile_height, channels, slices, frames);
	W = round(1.1*ncolumns*tile_width-(ncolumns-1)*overlap_x);
	H = round(1.1*nrows*tile_height-(nrows-1)*overlap_y);	
	newImage("Tiled", "32-bit grayscale-mode", W, H, channels, slices, 1);
	dst = getImageID();
	newImage("acc", "32-bit grayscale-mode", W, H, channels, slices, 1);
	acc = getImageID();
	newImage("tmp1", "32-bit grayscale-mode",  W, H, 1,1,1);
	tmp1 = getImageID();
	newImage("tmp2", "32-bit grayscale-mode",  W, H, 1,1,1);
	tmp2 = getImageID();
	b = 2; // border
	for (j = 0; j < nrows; j++) {
		for (i = 0; i < ncolumns; i++) {
			k = i + nrows * j;
			x = getResult("x", k);
			y = getResult("y", k);			
			for (channel = 1; channel <= channels; channel++) {
				for (slice = 1; slice <= slices; slice++) {	
					selectImage(src);
					Stack.setPosition(channel, slice, i+ncolumns*j+1);
					makeRectangle(b,b,tile_width-2*b,tile_height-2*b);
					run("Copy");
							
					// create a mask for the tile
					selectImage(tmp2);
					run("Select None"); setColor(0.0); fill();	
					makeRectangle(x+b,y+b,tile_width-2*b,tile_height-2*b);
					setColor(1.0); fill();											
					run("Macro...", "code=v=(v>0.0)");
					run("Minimum...", "radius="+overlap_x/4);
					run("Gaussian Blur...", "sigma="+overlap_x/8);				
						
					// add the tmp image to the accumulator acc
					selectImage(acc); Stack.setPosition(channel, slice, 1);
					imageCalculator("Add", acc, tmp2);
							
					// paste the tile in tmp1
					selectImage(tmp1);
					run("Select None");
					setColor(0); fill();
					makeRectangle(x+b,y+b,tile_width-2*b,tile_height-2*b);
					run("Paste");
					run("Select None");
					imageCalculator("Multiply", tmp1, tmp2);
					// add the tmp image to the accumulator dst
					selectImage(dst);Stack.setPosition(channel, slice, 1);
					imageCalculator("Add", dst, tmp1);
				}
			}
		}
	}	
	imageCalculator("Divide stack", dst, acc);
	selectImage(acc); close();
	selectImage(tmp1); close();
	selectImage(tmp2); close();
	selectImage(src);run("Select None");
	selectImage(dst);
	run("Select None");
	if (nSlices > 1) {
		Stack.setDisplayMode("composite");
		Stack.setActiveChannels("111");
		for (c = 1; c <= channels; c++){ Stack.setChannel(c); resetMinAndMax(); }
	}
	
	return dst;
}

function add_sprite(src,dst,cd,zd,td) {
	// not used
	selectImage(src);
	run("Select None");
	w = getWidth();
	h = getHeight();
	buf = newArray(w*h);
	for (y = 0; y < h; y++) {
		for (x = 0; x < w; x++) {
			buf[x+y*w] = getPixel(x, y);	
		}
	}
	selectImage(dst);
	Stack.setPosition(cd, zd, td);
	run("Select None");					
	for (y = 0; y < h; y++) {
		for (x = 0; x < w; x++) {
			v = getPixel(x,y);
			setPixel(x,y,v+buf[x+y*w]);	
		}
	}
}

function apply_transform(ncolumns,nrows,overlap_x,overlap_y) {
	// apply the shift in the result table to the stack
	// apply to all channels & slices
	src = getImageID();
	Stack.getDimensions(tile_width, tile_height, channels, slices, frames);
	for (k = 0; k < frames; k++) {
		dx0 = getResult("dx", k);
		dy0 = getResult("dy", k);
		for (c = 1; c <= channels; c++) {
			for (z = 1; z <= slices; z++) {
				Stack.setPosition(c, z, k+1);
				run("Translate...", "x="+dx0+" y="+dy0+" interpolation=Bicubic slice");
			}
		}
	}
}

function listOverlapingTiles(t0,t1,a0) {
	// list all overlaping tiles
	n = nResults;	
	c0 = getTileCoordsFromTable(t0);
	list = newArray(n);
	k = 0;
	for (i = 0; i < n; i++) if (i != t0 && i != t1) {
		ci = getTileCoordsFromTable(i);
		if (areaTilesOverlap(c0, ci) > a0) {
			list[k] = i;
			k++;
		}
	}	
	return Array.trim(list,k);
	
}

function align_tiles(ncolumns,nrows,overlap_x,overlap_y) {
	n = nResults;	
	j0 = log(overlap_x) / log(2) - 3;		
	print("min scale " + pow(2.0,j0));	
	for (iter = 0; iter < 10; iter++) {
		// select a random starting point
		t0 = Math.floor(random * n);			
		t1 = t0;
		t2 = t0;				
		for (j = j0; j >=0; j--) {
			delta = 0;
			nsteps = 0;
			for (step = 0; step < nrows*ncolumns; step++) {
				list = listOverlapingTiles(t0,t1,overlap_x*overlap_y);			
				if (list.length==0) {print("empty list");break;}
				t2 = t1;
				t1 = t0;			
				t0 = list[Math.floor(list.length*random)];								
				shift = ncropAndEstimateShift(t0,t1,pow(2.0,-j));	
				Stack.setFrame(t1+1);
				run("Select None");
				run("Translate...", "x="+-shift[0]+" y="+-shift[1]+" interpolation=Bicubic slice");
				dx0 = getResult("dx", t1);
				dy0 = getResult("dy", t1);			
				setResult("dx", t1, dx0-shift[0]);
				setResult("dy", t1, dy0-shift[1]);
				updateResults();
				delta += sqrt(shift[0]*shift[0]+shift[1]*shift[1]);
				nsteps++;
			}
			//print(delta/nsteps);
		}
	}
	run("Select None");
}

function getOverlapRectOnTile(id,i,j) {	
	// Duplicate overlapping region from tile at position i based on result table	
	A = getTileCoordsFromTable(i);	
	B = getTileCoordsFromTable(j);
	x0 = maxOf(A[0],B[0]);
	y0 = maxOf(A[1],B[1]);
	x1 = minOf(A[0]+A[2],B[0]+B[2]);
	y1 = minOf(A[1]+A[3],B[1]+B[3]);
	selectImage(id);	
	run("Select None");
	Stack.setFrame(i+1);
	makeRectangle(x0-A[0], y0-A[1], x1-x0, y1-y0);		
	run("Duplicate...", "title=I");
	run("Unsharp Mask...", "radius=2 mask=0.9");
	I = getImageID;
	selectImage(id);
	run("Select None");
	Stack.setFrame(j+1);
	makeRectangle(x0-B[0], y0-B[1], x1-x0, y1-y0);	
	run("Duplicate...", "title=J");
	run("Unsharp Mask...", "radius=2 mask=0.9");
	J = getImageID;	
	return newArray(I,J);
}


function ncropAndEstimateShift(t0,t1,scale) {
	id0=getImageID();
	run("Select None");
	id = getOverlapRectOnTile(id0,t0,t1);			
	shift = estimateShift(id[0],id[1],scale);
	selectImage(id[0]); close();
	selectImage(id[1]); close();
	if (isNaN(shift[0])) {shift[0]=0;}
	if (isNaN(shift[1])) {shift[1]=0;}
	return shift;
}


function estimateShift(id1,id2,scale) {
	sx = 0;
	sy = 0;
	n = 1;
	for (iter = 0; iter < 10 && n > 0.01; iter++) {
		selectImage(id2); run("Duplicate...", " "); id3 = getImageID();
		if (iter>0) {
			run("Translate...", "x="+-sx+" y="+-sy+" interpolation=Bilinear");
		}
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
	w = round(getWidth()*scale);
	h = round(getHeight()*scale);
	run("Duplicate...", " ");
	run("32-bit");
	if (scale != 1) {
		run("Size...", "width="+w+" height="+h+" depth=1 constrain average interpolation=Bilinear");	
	}
	ix = getImageID();rename("ix");
	selectImage(ix);
	run("Duplicate...", " ");
	iy = getImageID();rename("iy");

	imageCalculator("Subtract create 32-bit", id2, id1);
	it = getImageID();rename("it");
	if (scale > 0) {
		run("Size...", "width="+w+" height="+h+" depth=1 constrain average interpolation=Bilinear");	
	}

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
	return newArray(vx/scale,vy/scale);
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

function split_tiles(ncolumns,nrows,overlap_x,overlap_y,noise) {
	src = getImageID();
	Overlay.remove();
	Stack.getDimensions(width, height, channels, slices, frames);
	tile_width = Math.floor((width+(ncolumns-1)*overlap_x) / ncolumns);
	tile_height = Math.floor((height+(nrows-1)*overlap_y) / nrows);
	newImage("Tiles", "8-bit", tile_width, tile_height, channels, slices, ncolumns*nrows);
	dst = getImageID();
	k = 0;
	for (j = 0; j < nrows; j++) {
		for (i = 0; i < ncolumns; i++) {
			x = round(i * (tile_width - overlap_x) + noise*random);
			y = round(j * (tile_height - overlap_y) + noise*random);
			setResult("ROI", k,"Tile-"+(k+1));
			setResult("Position",k,k+1);
			setResult("x",k,x);
			setResult("y",k,y);			
			setResult("width",k,tile_width);
			setResult("height",k,tile_width);
			setResult("dx",k,0);
			setResult("dy",k,0);
			k++;
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
	selectImage(src);
	run("Select None");
	selectImage(dst);
	run("Select None");
	return dst;
}

//// OLDER CODE

function align_tiles_old(ncolumns,nrows,overlap_x,overlap_y) {
	src = getImageID();
	Stack.getDimensions(tile_width, tile_height, channels, slices, frames);
	ref = newArray(nrows*ncolumns);
	idx = newArray(nrows*ncolumns);
	sideref = newArray(nrows*ncolumns);
	sideidx = newArray(nrows*ncolumns);	
	for (k = 0; k < frames; k++) {
		setResult("x", k, 0);
		setResult("y", k, 0);
	}	
	createPath(ref,idx,sideref,sideidx,ncolumns,nrows,0);
	processPath(ref,idx,sideref,sideidx,overlap_x,overlap_y);
	createPath(ref,idx,sideref,sideidx,ncolumns,nrows,1);
	processPath(ref,idx,sideref,sideidx,overlap_x,overlap_y);
	createPath(ref,idx,sideref,sideidx,ncolumns,nrows,2);
	processPath(ref,idx,sideref,sideidx,overlap_x,overlap_y);
	createPath(ref,idx,sideref,sideidx,ncolumns,nrows,3);
	processPath(ref,idx,sideref,sideidx,overlap_x,overlap_y);
	createPath(ref,idx,sideref,sideidx,ncolumns,nrows,0);
	processPath(ref,idx,sideref,sideidx,overlap_x,overlap_y);
	run("Select None");
}

function processPath(ref,idx,sideref,sideidx,overlap_x,overlap_y) {
	// register images following the path defined by the arguments
	for (k = 1; k < ref.length; k++) {
		shift = cropAndEstimateShift(ref[k],idx[k],sideref[k],sideidx[k],overlap_x,overlap_y);
		Stack.setFrame(idx[k]);
		run("Select None");
		run("Translate...", "x="+-shift[0]+" y="+-shift[1]+" interpolation=Bicubic slice");
		dx0 = getResult("x", idx[k]-1);
		dy0 = getResult("y", idx[k]-1);
		setResult("x",idx[k]-1,dx0-shift[0]);
		setResult("y",idx[k]-1,dy0-shift[1]);
	}
}

function createPath(ref,idx,sideref,sideidx,ncolumns,nrows,type) {
	// create a path across tiles with matching sides
	k = 0;
	if (type==0) { // forward column major
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
	} else if(type==1) { // forward row major
		for (i = 0; i < ncolumns; i++) {
			for (j = 0; j < nrows; j++) {	
				if (k != 0) {
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
	} else if (type==2) { // backward column major
		for (j = nrows-1; j >= 0; j--) {
			for (i = ncolumns-1; i >= 0; i--) {
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
	} else if(type==3) { // backward row major
		for (i = ncolumns-1; i >= 0; i--) {
			for (j = nrows-1; j > = 0; j--) {	
				if (k != 0) {
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
	}
}

function duplicate_side(name,n,side,overlap_x,overlap_y) {
	// duplicate a band of the tile from one of the side 1:top, 2: right, 2: bottom, 3: left
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


function correct_illumination() {	
	// correct illumination 
	setBatchMode(true);
	run("Duplicate...", "title=y duplicate");
	run("32-bit");
	y = getImageID();
	run("Min...", "value=0.01 stack");
		
	selectImage(y);
	run("Z Project...", "projection=Median");
	run("32-bit");
	rename("m");
	run("Gaussian Blur...", "sigma=10");
	getStatistics(voxelCount, mean, min, max, stdDev);		
	run("Macro...", "code=v=0.25+0.75*(v-"+min+")/("+max+"-"+min+") stack");
	
	
	m = getImageID();

	selectImage(m);
	run("Duplicate...", "title=b");
	run("32-bit");
	setColor(0.0);
	fill();	
	b = getImageID();	

	imageCalculator("Subtract create 32-bit stack", y, b); rename("x"); x = getImageID(); 
	imageCalculator("Divide stack", x, m);	
	run("Min...", "value=0.01 stack");
	
	for (iter = 0; iter < 10; iter++) {				
		imageCalculator("Multiply create 32-bit stack", x, m); rename("tmp"); tmp=getImageID(); // mx
		imageCalculator("Add stack", tmp, b); // mx+b	
		imageCalculator("Divide stack", tmp, y); // (mx+b) / y
		selectImage(tmp);		
		run("Min...", "value=0.01 stack");
		run("Reciprocal","stack"); // y / (mx+b)		
		run("Z Project...", "projection=[Average Intensity]");rename("avg");avg=getImageID();
		imageCalculator("Multiply", m, avg); // update m = m*y/mx+b		
		selectImage(tmp); close();
		selectImage(avg); close();
		selectImage(m);run("Gaussian Blur...", "sigma=1");
		run("Min...", "value=0.01 stack");
		
		imageCalculator("Multiply create 32-bit stack", x, m); tmp = getImageID(); // mx
		imageCalculator("Add stack", tmp, b); // mx+b		
		imageCalculator("Divide stack", tmp, y); // (mx+b) / y
		selectImage(tmp);
		run("Min...", "value=0.01 stack");
		run("Reciprocal","stack"); // y / (mx+b)
		imageCalculator("Multiply", x, tmp);
		selectImage(tmp); close();		
	}	
	selectImage(x);
	setBatchMode(false);
}

function linear_regression_stack() {
	// linear regression along tiles y = alpha + beta x
	// beta = (n Sxy - SxSy) / (n Sxx-Sx^2) = (n Sxy/Sx - Sy) / ((n * Sxx-Sx^2)/Sx)
	// alpha = Sy/n - beta Sx / n = (Sy / Sx - beta) * Sx / n	
	setBatchMode(false);
	id0 = getImageID;
	run("Duplicate...", "duplicate");
	run("Gaussian Blur...", "sigma=10 stack");
	y = getImageID();
	// Compute Sx	
	Stack.getDimensions(width, height, channels, slices, n);
	Sx = n * (n + 1) / 2; // thank you euler
	Sxx = n * (n + 1) * ( 2 * n + 1) / 6;
	// compute Sy
	run("Z Project...", "projection=[Sum Slices]");run("32-bit");rename("Sy");sy=getImageID();
	
	// compute XY
	selectImage(y);
	run("Duplicate...", "title=XY duplicate");selectWindow("XY");run("32-bit");xy=getImageID();	
	Stack.getDimensions(width, height, channels, slices, n);
	for(i = 1; i <= n; i++){Stack.setFrame(i); setColor(i);fill();}	
	imageCalculator("Multiply stack", xy, y);
	// compute SXY
	run("Z Project...", "projection=[Sum Slices]");rename("Sxy");sxy=getImageID();
	// compute beta
	run("Multiply...", "value="+n/Sx+" stack");	// n Sxy / Sx
	imageCalculator("Subtract stack", sxy, sy); // n Sxy / Sx - Sy
	run("Divide...", "value="+(n*Sxx-Sx*Sx)/Sx+" stack");	// n Sxy / Sx
	rename("beta");
	// compute alpha
	selectImage(sy);
	run("Divide...", "value="+Sx+" stack"); // Sy / Sx
	imageCalculator("Subtract", sy, sxy); // Sy / Sx - beta
	run("Multiply...", "value="+Sx/n+" stack"); // (Sy / Sx - beta) * Sx /n
	rename("alpha");

	selectImage(y);close();
	selectImage(xy);close();
}

