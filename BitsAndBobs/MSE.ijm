#@Dataset id1
#@Dataset id2
/*
 * Compute the MSE between two images
 * 
 * Jerome for Akaash 2023
 */
imageCalculator("Subtract create 32-bit", id1,id2);
selectWindow("Result of blobs.gif");
run("Square");
run("Select None");
mse = getValue("Mean");
close();
print(mse);