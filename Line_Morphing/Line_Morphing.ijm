/*
 * Line Moprhing
 * 
 * Morph 2 multi line ROI at different z position 
 * 
 * Usage: 
 * - Add 2 ROI to the ROI manager
 * - Press Run
 */
id = getImageID();
print("\\Clear");
// get the coordinates of the two lines
p0 = getROICoordinates(0);
Roi.getPosition(channel, z0, frame);
p1 = getROICoordinates(1);
Roi.getPosition(channel, z1, frame);
//print(mat2str(p0,"p0"));

// compute paiwie distances
D = pairwiseDistance(p0,p1);
print(mat2str(D,"D"));
//mat2img(D,"D");

// find the mapping
P = sinkhorn_uniform(D,10);
//mat2img(P,"P");

// displacement interpolation from one line into the other
selectImage(id);
n = matcol(p0);
m = matcol(p1);
for (z = z0; z <= z1; z++) {
	t = (z-z0) / (z1-z0);
	x = newArray(n);
	y = newArray(n);
	for (i = 0; i < n; i++) {
		for (j = 0; j < m; j++) {
			if (matget(P,i,j) > 1e-3) {
				x[i] = (1-t) * matget(p0,0,i) + t * matget(p1,0,i);
				y[i] = (1-t) * matget(p0,1,i) + t * matget(p1,1,i);		
			}
		}
	}
	roiManager("select",0);
	Roi.setPolylineSplineAnchors(x, y);
	Stack.setSlice(z);
	roiManager("add");
	wait(100);
}


function sinkhorn_uniform(C,gamma) {
	// Optimal transport using uniform weights
	n = matrow(C);
	m = matcol(C);
	p = matones(n,1);
	q = matones(m,1);
	P = sinkhorn(C,p,q,gamma);
	return P;
}

function sinkhorn(C,p,q,gamma) {
	// Optimal transport 
	// N = [300,200]
	// gamma = 0.01
	// xi = exp(-C/gamma)
	// b = ones(N(2),1)
	// a = p ./ (xi*b)
	// b = q ./ (x'*a)
	// final is diag(a)*xi*diag(b)
	n = matrow(C);
	m = matcol(C);
	assert((matsize(p) == n), "Size of C ["+n+"x"+m+"] and p ["+matsize(p)+"] do not match");
	assert((matsize(q) == m), "Size of C ["+n+"x"+m+"] and p ["+matsize(q)+"] do not match");
	p = matreshape(p,n,1);
	q = matreshape(q,m,1);
	// normalize p and q
	p = mat_normalize_sum(p);
	q = mat_normalize_sum(q);
	// compute a default gamma value from C
	if (gamma < 0) {
		minC = C[2];
		for (i = 3; i < C.length; i++) {minc = minOf(minc, C[i]);}
		gamma = minC / 100.0;
	}
	// compute xi = exp(-C/gamma)
	xi = matzeros(n,m);	
	for (i = 0; i < n; i++) {
		for (j = 0; j < m; j++) {
			cij = matget(C,i,j);
			val = exp(-cij/gamma);
			matset(xi,i,j,val);
		}
	}
	xit = mattranspose(xi);
	b = matones(m,1);
	for (iter = 0; iter < 10; iter++) {
		xib = matmul(xi,b);
		a = mat_div(p,xib);
		xita = matmul(xit,a);
		b = mat_div(q,xita);
	}
	a = matdiag(a);
	b = matdiag(b);
	P = matmul(xi, b);
	P = matmul(a, P);
	return P;
}

function mat_normalize_sum(A) {
	// normalize A by the sum over all indices
	B = A;
	s = 0;
	for (i = 2; i < A.length; i++) {
		 s += A[i];
	}
	for (i = 2; i < A.length; i++) {
		 B[i] = A[i] / s;
	}
	return B;
}

function assert(cond, msg) {
	// check that condition is true
	if (!(cond)) {
		exit(msg);
	}
}

function pairwiseDistance(p,q) {
	// Compute the pairwie distance between p and q
	// P : 2xn coordinate matrix
	// Q : 2xm coordinate matrix
	// return distance [n x m] matrix
	n = matcol(p);
	m = matcol(q);
	print("n="+n);
	print("m="+m);
	
	D = matrow(p);
	a = matzeros(n,m);
	for (i = 0; i < n; i++) {
		for (j = 0; j < m; j++) {
			s = 0;
			for (k = 0; k < D; k++) {
				d = matget(p,k,i) - matget(q,k,j);
				s += d*d;
			}
			matset(a,i,j,s);
		}
	}
	return a;
}


function getROICoordinates(roi_id) {
	// get the spline anchor of the ROI
	// return a 2xn matrixs
	roiManager("select",roi_id);	
	Roi.getSplineAnchors(x, y);
	n = x.length;
	vals = Array.concat(x,y);
	m = matnew(2,n,vals);
	return m;
}


/////////////////////////////////////////////////////////
// Basic Linear Algebra
////////////////////////////////////////////////////////
function matnew(n,m,vals) {
	// Initialize a matrix with values
	// Column major, N rows and M columns
	// Matrices are array with 2 first elements being the dimensions
	A = newArray(2+n*m);
	A[0] = n;
	A[1] = m;
	for (i = 0; i < vals.length; i++) {
		A[2+i] = vals[i];
	}
	return A;
}

function matrow(A) {
	// return the number of rows
	return A[0];
}

function matcol(A) {
	// return the number of columns
	return A[1];
}

function matsize(A) {
	// return the number of elements of A
	return A.length - 2;
}

function matzeros(n,m) {
	// Return an NxM matrix initialized with zeros
	A = newArray(2+n*m);
	A[0] = n;
	A[1] = m;
	return A;
}

function matones(n,m) {
	// Return an NxM matrix initialized with zeros
	A = matzeros(n,m);
	for (i = 2 ; i < A.length; i++) {
		A[i] = 1;
	}
	return A;
}

function matreshape(A,n,m) {
	// Reshape the matrix
	B = A;
	B[0] = n;
	B[1] = m;
	return B;
}

function matdiag(v) {
	// create a diagonal matrix with the values of the vector v
	n = v.length - 2;
	A  = matzeros(n,n);
	for (i = 0 ; i < n; i++) {
		matset(A,i,i,matget(v,i,0));
	}
	return A;
}

function matset(A,i,j,val) {
	// set value of A[i,j] in place
	if (i >= 0 && i < A[0] && j >=0 && j < A[1]) {
		A[2+j+A[1]*i] = val;
	} else {
		print("Error in matset: out of bound [" + i + ","+j+"] in A["+A[0]+","+A[1]+"]");
		please_check_the_message_in_the_log();
	}
}

function matget(A,i,j) {
	// get the value A[i,j]
	if (i >= 0 && i < A[0] && j >=0 && j < A[1]) {
		return A[2+j+A[1]*i];
	} else {
		print("Error in matset: out of bound [" + i + ","+j+"] in A["+A[0]+","+A[1]+"]");
		please_check_the_message_in_the_log();
	}
}

function mattranspose(A) {
	// transpose matrix A
	M = matcol(A);
	N = matrow(A);
	B = matzeros(M,N);
	for (n=0;n<N;n++) {
		for (m=0;m<M;m++) {
			B[2+n+m*N] = A[2+m+n*M];
		}
	}
	return B;
}

function mat2str(A,name) {
	// Convert the matrix to a string
	if (A.length == 0){
		str = name + "is an empty matrix";
	} else {
		N = A[0];
		M = A[1];
		str = name + " ["+N+"x"+M+"] :\n";
		for (n=0;n<N;n++) {
			for (m=0;m<M;m++) {
				str = str + A[2+m+n*M];
				if ( m == M-1 ){ str = str + "\n";}
				else  {str = str + ", ";}
			}
		}
	}
	return str;
}

function mat2img(A,name) {
	// create an image and fill it with the matrix values
	width = matcol(A);
	height = matrow(A);
	newImage(name, "32-bit", width, height, 1);
	id = getImageID();
	for (y = 0; y < height; y++) {
		for (x = 0; x < width; x++) {
			setPixel(x, y, matget(A,y,x));
		}
	}
	run("Enhance Contrast", "saturated=0.0");
	run("In [+]");
	run("In [+]");
	return id;
}

function matadd(A,B,alpha) {
	// Add 2 matrices A and B (A + alpha B)
	C = matzeros(A[0],A[1]);
	for (i = 2; i <A.length; i++){C[i] = A[i]+alpha*B[i];}
	return C;
}

function matrep(A,K,L) {
	// repeat matrix A
	N = A[0];
	M = A[1];
	B = matzeros(N*K,L*M);
	for (k = 0; k < K; k++) {
		for (l = 0; l < L; l++) {	
			for (n = 0; n < N; n++) {
				for (m = 0; m < M; m++) {
					B[2+m+l*M+(M*L)*(n+k*N)] = A[2+m+n*M];
				}
			}
		}
	}
	return B;
}

function matmul(A,B) {
	// matrix multiplication C[nxk] = A[nxm] x B[mxk] 
	// C[n,k] = sum A[n,m] B[m,k]
	// (column major matrices)
	N = A[0];
	M = A[1];
	K = B[1];
	if (B[0] != A[1]) {
		print("matmul(A,B): incompatible dimensions A ["+ A[0]+"x"+A[1]+"], B ["+B[0]+"x"+B[1]+"].");
		please_check_the_message_in_the_log();
		
	}
	C = matzeros(N,K);
	for (n = 0; n < N; n++) {
		for (k = 0; k < K; k++) {
			for (m = 0; m < M; m++) {
				C[2+k+n*K] = C[2+k+n*K] + A[2+m+n*M] * B[2+k+m*K];
			}
		}
	}
	return C;
}

function mat_mul(A,B) {
	// multiply entry wise the two matrices
	// C[i,j] = A[i,j] * B[i,j]
	assert((A.length == B.length), "A and B have incompatible size");
	N = A[0];
	M = A[1];
	C = matzeros(N,M)
	for (i = 2; i < C.length; i++) {
		C[i] = A[i] * B[i];
	}
	return C;
}

function mat_div(A,B) {
	// divide entry wise A by B
	// C[i,j] = A[i,j] / B[i,j]
	assert((A.length == B.length), "A and B have incompatible size");
	N = A[0];
	M = A[1];
	C = matzeros(N,M);
	for (i = 2; i < C.length; i++) {
		C[i] = A[i] / B[i];
	}
	return C;
}


function matsum(A,dim) {
	// sum values of A algon 1 dimension (0 or 1)
	N = A[0];
	M = A[1];
	if(dim==0) {
		S = matzeros(1,M);
		for (n = 0; n < N; n++) {
			for (m = 0; m < M; m++) {
				S[2+m] = S[2+m] + A[2+m+n*M];
			}
		}
	} else {
		S = matzeros(N,1);
		for (n = 0; n < N; n++) {
			for (m = 0; m < M; m++) {
				S[2+n] = S[2+n] + A[2+m+n*M];
			}
		}
	}
	return S;
}

function matmean(A,dim) {
	// average along one dimension
	B = matsum(A,dim);
	for (i = 2; i < B.length; i++) {
		B[i] = B[i] / A[dim];
	}
	return B;
}

function solve(A,B) {
	// Solve A*X=B by iterating X=X-At(AX - B)
	// A[NxM] x X[MxK] = B [NxK]
	print(mat2str(A,"A"));
	print(mat2str(B,"B"));
	N = A[0];
	M = A[1];
	K = B[1];
	At = mattranspose(A);
	AtA = matmul(At,A);
	AtB = matmul(At,B);
	X = matzeros(M,K);
	a = matnorm(AtA);
	for (iter = 0; iter < 1000; iter++) {
		AtAX = matmul(AtA,X);
		dX = matadd(AtAX,AtB,-1);
		X = matadd(X, dX, -1/(2*a));
	}
	return X;
}

function matnorm(A) {
	// L2 norm of the matrix (sum a_ij^2)
	n = 0;
	for (i = 2; i < A.length; i++) {
		n = n + A[i]*A[i];
	}
	return sqrt(n);
}

