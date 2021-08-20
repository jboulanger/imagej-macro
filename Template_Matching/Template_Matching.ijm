
// Load the pixels of the image in an array
I = convertImageToArray();
// Load the pixels of the template in an array
run("Duplicate...", " ");
T = convertImageToArray();
close();
C = correlateTemplate(I,T);
createImageFromArray("correlation","32-bit",C);

function matchTemplateNormXCorr(I,J) {
	/* Correlate the image I J with template J
	 * Return the normalized cross correlation map
	 */
	wi = I[0];
	hi = I[1];
	wj = T[0];
	hj = T[1];
	C = newArray(2+wi*hi);
	C[0] = wi;
	C[1] = hi;
	mui = computeImageMean(I);
	mut = computeImageMean(T);
	for (yi = hj; yi < hi-hj; yi++) {
		showProgress(yi-hj,hi-2*hj);
		for (xi = wj; xi < wi-wj; xi++) {
			s1 = 0;
			s2 = 0;
			s3 = 0;
			for (yj = 0; yj < hj; yj++) {
				for (xj = 0; xj < wj; xj++) {				
					vi = I[2+xi+xj+(wi)*(yi+yj)] - mui;
					vj = J[2+xj+(wj)*yj] - mut;
					s1 += vi * vj;					
					s2 += vi * vi;
					s3 += vj * vj;
				}	
			}
			C[2+xi+wi*yi] = s1;
		}
	}
	return C;
}

function matchTemplateSAD(I,J) {
	/* Correlate the image I J with template J
	 * Return the normalized cross correlation map
	 */
	wi = I[0];
	hi = I[1];
	wj = T[0];
	hj = T[1];
	C = newArray(2+wi*hi);
	C[0] = wi;
	C[1] = hi;	
	for (yi = hj; yi < hi-hj; yi++) {
		showProgress(yi-hj,hi-2*hj);
		for (xi = wj; xi < wi-wj; xi++) {
			s1 = 0;			
			for (yj = 0; yj < hj; yj++) {
				for (xj = 0; xj < wj; xj++) {				
					vi = I[2+xi+xj+(wi)*(yi+yj)];
					vj = J[2+xj+(wj)*yj];
					s += abs(vi - vj);
				}	
			}
			C[2+xi+wi*yi] = s / (wj*hj);
		}
	}
	return C;
}

function computeImageMean(A){
	s = 0;
	for (i=2;i<A.length;i++) {s+=A[i];}
	return s / (A.length-2);
}


function createImageFromArray(title,type,A) {
	/* Create a 2D image from the array A
	 * 
	 */
	width = A[0];
	height = A[1];
	newImage(title, type, width, height, 1);
	for (y = 0; y < height; y++) {
		for (x = 0; x < width; x++) {
			setPixel(x,y,A[2+x + width * y]);
		}
	}
}


function convertImageToArray() {
	/* Convert the current image to an array
	 * 
	 */
	width = getWidth();
	height = getHeight();
	A = newArray(2+width*height);
	A[0] = width;
	A[1] = height;
	for (y = 0; y < height; y++) {
		for (x = 0; x < width; x++) {
			A[2+x + width * y] = getPixel(x, y);
		}
	}
	return A;
}

