/*
 * segment and compute overlap of masks
*/
run("Duplicate...", "duplicate");
run("32-bit");
run("Square Root", "stack");
run("Median...", "radius=1 stack");
run("Subtract Background...", "rolling=10 stack");
id = getImageID();

run("Duplicate...", "duplicate channels=1");
id1 = getImageID();
Stack.getStatistics(voxelCount, mean, min, max, stdDev);
setThreshold(mean + 4 * stdDev, max);
run("Make Binary", "background=Dark");
t1 = getTitle();

selectImage(id);
run("Duplicate...", "duplicate channels=2");
id2 = getImageID();
Stack.getStatistics(voxelCount, mean, min, max, stdDev);
setThreshold(mean + 4 * stdDev, max);
run("Make Binary", "background=Dark");
t2 = getTitle();

imageCalculator("And create stack", id1, id2);
t3 = getTitle();

run("Merge Channels...", "c1=["+t1+"] c2=["+t2+"] c3=["+t3+"] create");

