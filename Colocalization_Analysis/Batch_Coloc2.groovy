/**
 * Coloc_Analysis.groovy
 *
 * Performs pixel-based colocalization analysis between two channels of a
 * multichannel image. Computes Pearson's correlation coefficient, Manders'
 * overlap coefficients.
 *
 * Usage:
 *   Run via Plugins > Macros > Run... in Fiji.
 *
 * Dependencies:
 *   - Fiji with ImgLib2 (net.imglib2)
 *
 * Author:  Jerome Boulanger
 * Date:    11 March 2026
 * Version: 1.
 */

// -- Script Parameters -------------------------------------------------------

#@ File    (label="Select Multichannel Image")          imageFile
#@ Integer (label="Channel A", value=1)                 chIndexA
#@ Integer (label="Channel B", value=2)                 chIndexB
#@ Double  (label="Threshold A", value=0)               thresholdA
#@ Double  (label="Threshold B", value=0)               thresholdB
#@ DatasetIOService io
#@ UIService        ui

// -- Imports -----------------------------------------------------------------

import net.imglib2.view.Views
import net.imglib2.RandomAccessibleInterval
import ij.measure.ResultsTable
import net.imglib2.util.Util
import ij.IJ
import java.io.PrintWriter
import java.io.StringWriter
import net.imglib2.util.Intervals




// ----------------------------------------------------------------------------
// computeColocalization
//
// Single-pass Pearson, Manders M1/M2, k1/k2 and MOC via Welford's algorithm.
// ----------------------------------------------------------------------------
Map computeColocalization(RandomAccessibleInterval img1, RandomAccessibleInterval img2,
                          double threshold1 = 0.0, double threshold2 = 0.0) {

    def c1 = Views.iterable(img1).cursor()
    def c2 = Views.iterable(img2).cursor()

    long   n      = 0
    double mean1  = 0, mean2  = 0
    double cov    = 0
    double var1   = 0, var2   = 0
    double sum1   = 0, sum2   = 0
    double col1   = 0, col2   = 0
    double dot12  = 0
    double sqSum1 = 0, sqSum2 = 0

    while (c1.hasNext()) {
        double v1 = c1.next().getRealDouble()
        double v2 = c2.next().getRealDouble()
        n++
        double d1 = v1 - mean1
        double d2 = v2 - mean2
        mean1 += d1 / n;   mean2 += d2 / n
        cov   += d1 * (v2 - mean2)
        var1  += d1 * (v1 - mean1);  var2 += d2 * (v2 - mean2)
        sum1  += v1;  sum2  += v2
        dot12 += v1 * v2
        sqSum1 += v1 * v1;  sqSum2 += v2 * v2
        if (v2 > threshold2) col1 += v1
        if (v1 > threshold1) col2 += v2
    }

    return [
        pearson : cov  / Math.sqrt(var1 * var2),
        m1      : col1 / sum1,
        m2      : col2 / sum2,
        k1      : dot12 / sqSum1,
        k2      : dot12 / sqSum2,
        moc     : dot12 / Math.sqrt(sqSum1 * sqSum2),
    ]
}


// -- Main Script Body --------------------------------------------------------

def imagePath = imageFile.getAbsolutePath()
println "[Coloc] Processing: ${imagePath}"
println "[Coloc] Channel A: ${chIndexA}  Channel B: ${chIndexB}"

def image = io.open(imagePath)
def imgA  = Views.hyperSlice(image.getImgPlus(), 2, chIndexA - 1)
def imgB  = Views.hyperSlice(image.getImgPlus(), 2, chIndexB - 1)


try {
    def coloc = computeColocalization(imgA, imgB, thresholdA, thresholdB)


    def rt = ResultsTable.getResultsTable() ?: new ResultsTable()
    rt.incrementCounter()
    rt.addValue("Image",                         imageFile.getName())
    rt.addValue("chA",                       	 chIndexA)
    rt.addValue("chB",                       	 chIndexB)
    rt.addValue("Threshold A",                   thresholdA)
    rt.addValue("Threshold B",                   thresholdB)
    rt.addValue("Pearson R",                     coloc.pearson)
    rt.addValue("Manders k1",                    coloc.k1)
    rt.addValue("Manders k2",                    coloc.k2)
    rt.addValue("Manders M1",                    coloc.m1)
    rt.addValue("Manders M2",                    coloc.m2)
    rt.addValue("Manders MOC",                   coloc.moc)
    rt.show("Coloc Results")

    println "[Coloc] Done."

} catch (Exception e) {
    def sw = new StringWriter()
    e.printStackTrace(new PrintWriter(sw))
    IJ.log("[Coloc] ERROR:\n${sw.toString()}")
}
