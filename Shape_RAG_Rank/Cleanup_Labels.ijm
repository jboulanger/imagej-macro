#@File("Input label image",style="open") inputFile
/*
 * Cleanup label image
 * 
 * - Remove border labels
 * - Remap labels
 * 
 * 
 * Requires MorpholibJ
 * 
 * Jerome Boulanger for Nanami and Vivi
 */

setBatchMode("hide");
open(inputFile);
run("Remove Border Labels", "left right top bottom");
run("Remap Labels");
ofile = File.getDirectory(inputFile) 
ofile += File.separator() + "LabelEdited" + File.separator()
ofile += File.nameWithoutExtension(inputFile)
ofile += "_edited.tiff";
print(ofile);
saveAs("TIFF", ofile);
setBatchMode("exit and display");
