%% Load the image and display the result
A = imread('test.tif',1);
B = imread('test.tif',2);
R = imref2d(size(A),[-1,1],[-1,1]);
imshowpair(A,B);
%% Load the deformation and apply to the second image
filename = '~/Desktop/Models.csv';
tab = readtable(filename);
switch tab.M12(2)
    case 3        
        M = reshape(tab.M12(3:end)',3,2)';                 
        tform = affine2d(M(:,[2,3,1]));
    case 6                    
        M =  reshape(tab.M12(3:end)',6,2)';        
        tform = images.geotrans.PolynomialTransformation2D(M(1,[1,2,3,5,4,6]),M(2,[1,2,3,5,4,6]));        
    otherwise
        error('Not implemented')
end
[Bw,Rw] = imwarp(B,R,tform,'OutputView',R);
imshowpair(A,Bw);