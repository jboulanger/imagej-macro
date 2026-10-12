/**
 * ============================================================================
 * 3D MULTI-CHANNEL SPOT DETECTION AND SPATIAL COLOCALIZATION PIPELINE
 * ============================================================================
 *
 * Overview:
 * ---------
 * This script provides an automated ImageJ/Fiji workflow for detecting 3D punctate
 * signals (spots/blobs) across multi-channel fluorescence microscopy images and 
 * quantifying their spatial colocalization across combinatorial channel sets.
 *
 * Key Pipeline Stages:
 * 1. Image I/O & Masking: Handles input images, 2D/3D stacks, CSV spot inputs,
 *    or synthetic test dataset generation.
 * 2. Feature Enhancement: Applies Difference of Gaussians (DoG), Laplacian of
 *    Gaussians (LoG), or 3D Top-Hat filtering to enhance diffraction-limited spots.
 * 3. Local Maxima & Thresholding: Identifies candidate centers using 3D local
 *    maxima combined with adaptive (MAD-based) or absolute intensity thresholding.
 * 4. Watershed Segmentation & Volume Filtering: Uses MorphoLibJ marker-controlled 
 *    3D watershed to filter out structures exceeding user-defined volume caps.
 * 5. Subpixel Localization: Refines spot centroids in continuous physical space
 *    using 1D orthogonal Gaussian intensity profile fitting.
 * 6. Spatial Colocalization Engine: Employs ImgLib2 KD-Trees for fast 3D radius
 *    neighbor searches (d_max) and evaluates power-set channel overlap codes.
 *
 * Dependencies (Fiji Plugins):
 * ----------------------------
 * - MorphoLibJ (inra.ijpb.*)
 * - ImgLib2 (net.imglib2.*)
 * - Bio-Formats (loci.plugins.BF)
 *
 * Developer: Research Software Engineering Core
 * ============================================================================
 */

// ============================================================================
// SCRIPT PARAMETERS (ImageJ Scripting UI Annotations)
// ============================================================================

#@ File (label="Input", value="image", description="'image' = current image, 'test' = generated test image, or a path to an image / a .csv coordinates file") inputPath
#@ String (label="Channels", value="1,2,3", description="Comma-separated list of 1-based channel indices (e.g. '1,2,3')") channelsStr
#@ String (label="Spot Size", value="0.5,0.5,0.5", description="Comma-separated expected spot radii in microns along XY per channel") spotSizeStr
#@ String (label="Max Spot Size", value="1,1,1", description="Comma-separated max allowed spot radii in microns; use -1 to disable volume filtering") maxSpotSizeStr
#@ String (label="Feature", choices={"DoG","LoG","Top hat"}, description="3D feature enhancement filter algorithm") feature
#@ String (label="Specificity[-log10]", value="3,3,3", description="Comma-separated specificity (-log10 PFA) per channel") specificityStr
#@ Boolean (label="Adaptive threshold", value=true, description="Use robust image statistics (MAD) for automatic thresholding") adaptive
#@ String (label="Mask channel", value="none", description="1-based channel index to use as a binary mask, or 'none'") maskStr
#@ Double (label="Proximity threshold [um]", value=0.5, description="Maximum Euclidean distance in physical units to consider two spots colocalized") dmax
#@ Boolean (label="Subpixel localization", value=true, description="Enable 3D continuous subpixel Gaussian refinement") subpixel
#@ Boolean (label="Close on exit", value=false, description="Automatically close image windows after workflow completion") closeOnExit
#@ Boolean (label="Save coordinates", value=false, description="Export localized spot coordinates to a CSV file") saveCoordinates
#@ Boolean (label="Z project", value=false, description="Generate a 2D Maximum Intensity Projection with overlay ROIs") doZProject



// Groovy language import
import groovy.transform.CompileStatic
import java.io.*
import java.util.stream.IntStream
import java.util.*

import java.awt.Color
import java.awt.Frame

// Standard Java & I/O Imports
import ij.*
import ij.gui.*
import ij.measure.*
import ij.plugin.*
import ij.plugin.filter.*
import ij.process.*
import ij.text.TextWindow

import inra.ijpb.binary.*
import inra.ijpb.distance.*
import inra.ijpb.label.*
import inra.ijpb.plugins.*
import inra.ijpb.watershed.*
import inra.ijpb.morphology.MinimaAndMaxima
import inra.ijpb.morphology.MinimaAndMaxima3D


import net.imglib2.algorithm.gauss3.Gauss3
import net.imglib2.algorithm.gradient.PartialDerivative
import net.imglib2.converter.Converters
import net.imglib2.Cursor
import net.imglib2.KDTree
import net.imglib2.RealPoint
import net.imglib2.RandomAccessibleInterval
import net.imglib2.type.numeric.real.FloatType
import net.imglib2.type.numeric.real.DoubleType
import net.imglib2.img.display.imagej.ImageJFunctions
import net.imglib2.interpolation.randomaccess.NLinearInterpolatorFactory
import net.imglib2.loops.LoopBuilder
import net.imglib2.neighborsearch.RadiusNeighborSearchOnKDTree
import net.imglib2.view.Views
import net.imglib2.algorithm.neighborhood.RectangleShape
import net.imglib2.algorithm.neighborhood.Neighborhood
import net.imglib2.img.display.imagej.ImageJFunctions

import loci.plugins.BF


// ============================================================================
// I/O & PARSING HELPER FUNCTIONS
// ============================================================================

/** Determines execution mode based on input string path. */
String getMode(String pathStr) {
    if (pathStr =~ /.*test/) return "test"
    if (pathStr =~ /.*image/) return "image"
    if (pathStr =~ /.*csv/) return "csv"
    return "file"
}

/** Parses comma-separated string to list of Integers. */
List parseCSVInt(String s) { 
    s.split(",")*.trim().findAll { it }*.toInteger() 
}

/** Parses comma-separated string to list of Doubles. */
List parseCSVFloat(String s) { 
    s.split(",")*.trim().findAll { it }*.toDouble() 
}

/** Parses optional binary mask channel parameter. */
Integer parseMaskStr(String s) { 
    (s != null && s.trim().isInteger()) ? s.trim().toInteger() : null 
}

/** Removes file extension from filename string. */
String stripExtension(String name) {
    int i = name.lastIndexOf(".")
    return i > 0 ? name.substring(0, i) : name
}

// ============================================================================
// VOLUME MANIPULATION & MATH HELPERS (ImgLib2 Accelerated)
// ============================================================================

/** Extracts 3D physical calibration voxel dimensions [dx, dy, dz] in microns. */
double[] voxelSize(ImagePlus imp) {
    Calibration cal = imp.getCalibration()
    return [cal.pixelWidth, cal.pixelHeight, cal.pixelDepth] as double[]
}

/** Compute the volume in microns */
double computePhysicalVolume(ImagePlus imp) {
    double[] v = voxelSize(imp) 
    
    int slices = imp.getNSlices() > 1 ? imp.getNSlices() : imp.getStackSize()
    double physicalArea = (imp.getWidth() * v[0]) * (imp.getHeight() * v[1])
    
    if (slices > 1) {
        return physicalArea * (slices * v[2]) // Volume in um^3
    }
    return physicalArea // Area in um^2
}


/** Extracts a single 1-based channel from an ImagePlus volume as 32-bit FloatProcessor stack. */
ImagePlus extractChannelFloat(ImagePlus imp, int channel) {
    ImageStack src = imp.getStack()
    ImageStack dst = new ImageStack(imp.getWidth(), imp.getHeight())
    for (int z = 1; z <= imp.getNSlices(); z++) {
        ImageProcessor ip = src.getProcessor(imp.getStackIndex(channel, z, 1))
        dst.addSlice(ip instanceof FloatProcessor ? ip.duplicate() : ip.convertToFloat())
    }
    ImagePlus out = new ImagePlus("channel" + channel, dst)
    out.setCalibration(imp.getCalibration().copy())
    return out
}

/** Deep duplicates a 32-bit float ImagePlus volume with calibration. */
ImagePlus duplicateFloat(ImagePlus imp, String title) {
    ImageStack s = imp.getStack()
    ImageStack d = new ImageStack(imp.getWidth(), imp.getHeight())
    for (int i = 1; i <= s.getSize(); i++) d.addSlice(s.getProcessor(i).duplicate())
    ImagePlus out = new ImagePlus(title, d)
    out.setCalibration(imp.getCalibration().copy())
    return out
}

/** Applies element-wise square root in-place. */
@CompileStatic
void sqrtInPlace(ImagePlus imp) {
    ImageStack stack = imp.getStack()
    int numSlices = stack.getSize()
    
    for (int i = 1; i <= numSlices; i++) {
        float[] pixels = (float[]) stack.getPixels(i)
        int len = pixels.length
        
        for (int k = 0; k < len; k++) {
            float v = pixels[k]
            pixels[k] = v > 0f ? (float) Math.sqrt(v) : 0f
        }
    }
}

/** Scales all voxel intensities by constant factor f in-place. */
@CompileStatic
void scaleAll(ImagePlus imp, float f) {
    ImageStack stack = imp.getStack()
    int numSlices = stack.getSize()
    
    for (int i = 1; i <= numSlices; i++) {
        float[] pixels = (float[]) stack.getPixels(i)
        int len = pixels.length
        
        for (int k = 0; k < len; k++) {
            pixels[k] *= f
        }
    }
}

/** Scales a specific slice's pixel values by constant factor f. */
@CompileStatic
void scaleSlice(ImagePlus imp, int slice, float f) {
    float[] px = (float[]) imp.getStack().getPixels(slice)
    for (int k = 0; k < px.length; k++) px[k] *= f
}

/** Element-wise multiplication: a = a * b */
@CompileStatic
void multiplyInPlace(ImagePlus a, ImagePlus b) {
    ImageStack sa = a.getStack()
    ImageStack sb = b.getStack()
    int numSlices = sa.getSize()
    
    for (int i = 1; i <= numSlices; i++) {
        float[] pxA = (float[]) sa.getPixels(i)
        float[] pxB = (float[]) sb.getPixels(i)
        int len = pxA.length
        
        for (int k = 0; k < len; k++) {
            pxA[k] *= pxB[k]
        }
    }
}

/** Element-wise subtraction: a = a - b */
@CompileStatic
void subtractInPlace(ImagePlus a, ImagePlus b) {
    ImageStack sa = a.getStack()
    ImageStack sb = b.getStack()
    int numSlices = sa.getSize()
    
    for (int i = 1; i <= numSlices; i++) {
        float[] pxA = (float[]) sa.getPixels(i)
        float[] pxB = (float[]) sb.getPixels(i)
        int len = pxA.length
        
        for (int k = 0; k < len; k++) {
            pxA[k] -= pxB[k]
        }
    }
}

/** Element-wise subtraction: a = a - b */
@CompileStatic
void subtractInPlaceImgLib2(
        RandomAccessibleInterval target, 
        RandomAccessibleInterval subtractor) {

    int splitDim = target.numDimensions() - 1
    int numSlices = (int) target.dimension(splitDim)

    if (numSlices > 4) {
        IntStream.range(0, numSlices).parallel().forEach { int i ->
            Cursor c1 = Views.iterable(Views.hyperSlice(target, splitDim, (long) i)).cursor()
            Cursor c2 = Views.iterable(Views.hyperSlice(subtractor, splitDim, (long) i)).cursor()
            while (c1.hasNext()) {
                FloatType v1 = (FloatType) c1.next()
                FloatType v2 = (FloatType) c2.next()
                v1.set(v1.get() - v2.get())
            }
        }
    } else {
        Cursor c1 = Views.iterable(target).cursor()
        Cursor c2 = Views.iterable(subtractor).cursor()
        while (c1.hasNext()) {
            FloatType v1 = (FloatType) c1.next()
            FloatType v2 = (FloatType) c2.next()
            v1.set(v1.get() - v2.get())
        }
    }
}

/** Binarizes volume in-place: pixel = 1.0 if pixel >= t else 0.0 */
@CompileStatic
void binarizeGE(ImagePlus imp, double t) {
    ImageStack s = imp.getStack()
    for (int i = 1; i <= s.getSize(); i++) {
        ImageProcessor ip = s.getProcessor(i)
        for (int k = 0; k < ip.getPixelCount(); k++) {
            ip.setf(k, (float) (ip.getf(k) >= t ? 1.0 : 0.0))
        }
    }
}

/** Inverts binary image: converts 0.0 -> 1.0 and non-zero -> 0.0 */
@CompileStatic
void zeroToOne(ImagePlus imp) {
    ImageStack s = imp.getStack()
    for (int i = 1; i <= s.getSize(); i++) {
        ImageProcessor ip = s.getProcessor(i)
        for (int k = 0; k < ip.getPixelCount(); k++) {
            ip.setf(k, (float) (ip.getf(k) == 0f ? 1.0 : 0.0))
        }
    }
}

/** Computes global minimum and maximum intensity values across the volume. */
@CompileStatic
double[] minMax(ImagePlus imp) {
    double mn = Double.MAX_VALUE, mx = -Double.MAX_VALUE
    ImageStack s = imp.getStack()
    for (int i = 1; i <= s.getSize(); i++) {
        float[] px = (float[]) s.getPixels(i)
        for (int k = 0; k < px.length; k++) {
            float v = px[k]
            if (v != v) continue
            if (v < mn) mn = v
            if (v > mx) mx = v
        }
    }
    return [mn, mx] as double[]
}

/** Performs 3D anisotropic Gaussian smoothing using ImgLib2 Gauss3. */
void gaussian(ImagePlus imp, double sx, double sy, double sz) {
    def img = ImageJFunctions.wrapReal(imp)
    double[] sigmas = (imp.getStackSize() > 1) ? ([sx, sy, sz] as double[]) : ([sx, sy] as double[])
    Gauss3.gauss(sigmas, Views.extendBorder(img), img)
}

/** Applies 3D Minimum or Maximum rank filter. */
ImagePlus rankFilter(ImagePlus imp, boolean max, double rx, double ry, double rz) {
    ImagePlus out
    if (imp.getStackSize() > 1) {
        ImageStack r = Filters3D.filter(imp.getStack(), max ? Filters3D.MAX : Filters3D.MIN,
            (float) rx, (float) ry, (float) rz)
        out = new ImagePlus(imp.getTitle(), r)
    } else {
        ImageProcessor ip = imp.getProcessor().duplicate()
        new RankFilters().rank(ip, rx, max ? RankFilters.MAX : RankFilters.MIN)
        out = new ImagePlus(imp.getTitle(), ip)
    }
    out.setCalibration(imp.getCalibration().copy())
    return out
}

/** Retrieves window reference by title or returns active window. */
ImagePlus resultImage(String title) {
    ImagePlus r = WindowManager.getImage(title)
    return r != null ? r : IJ.getImage()
}

/** Closes ImagePlus reference safely without prompting save dialogs. */
void dispose(ImagePlus imp) {
    if (imp == null) return
    imp.changes = false
    imp.close()
}

/** Close all image windows **/
void closeAllImages() {
    int[] ids = WindowManager.getIDList()
    if (ids != null) {
        for (int id : ids) {
            def imp = WindowManager.getImage(id)
            if (imp != null) {
                imp.changes = false
                imp.close()
            }
        }
    }
}

// ============================================================================
// STATISTICAL ESTIMATORS & AUTOMATIC THRESHOLDING
// ============================================================================

/** Normal Cumulative Distribution Function approximation. */
double normcdf(double x) {
    double s = 0, n = 1
    for (int i = 0; i < 100; i++) {
        n *= (2 * i + 1)
        s += Math.pow(x, 2 * i + 1) / n
    }
    return 0.5 + 1.0 / Math.sqrt(2.0 * Math.PI) * Math.exp(-0.5 * x * x) * s
}

/** Inverse Normal CDF (Percent Point Function) via Newton-Raphson. */
double normppf(double p) {
    if (Math.log10(p) < -10) {
        return -Math.sqrt(-4.2 * Math.log10(p) + 4.5)
    }
    double x = 0
    for (int i = 0; i < 50; i++) {
        double delta = (normcdf(x) - p) / Math.exp(-0.5 * x * x) * Math.sqrt(2 * Math.PI)
        x = x - delta
        if (Math.abs(delta) < 1e-6) break
    }
    return x
}

/** Computes robust statistics: Median and Median Absolute Deviation (MAD). */
double[] robustStats(ImagePlus imp, double mn, double mx) {
    int nbins = 1000
    if (mx <= mn) return [mn, 0d] as double[]
    ImageStack s = imp.getStack()

    long[] hist = new long[nbins]
    long n = 0
    double scale = nbins / (mx - mn)
    for (int i = 1; i <= s.getSize(); i++) {
        float[] px = (float[]) s.getPixels(i)
        for (int k = 0; k < px.length; k++) {
            float v = px[k]
            if (v != v) continue
            int b = (int) ((v - mn) * scale)
            if (b >= nbins) b = nbins - 1
            hist[b]++
            n++
        }
    }
    long acc = 0
    int mb = 0
    for (int i = 0; i < nbins; i++) {
        acc += hist[i]
        if (acc > n / 2) { mb = i; break }
    }
    double med = mn + (mx - mn) * mb / nbins

    double maxDev = Math.max(Math.abs(mn - med), Math.abs(mx - med))
    if (maxDev == 0) return [med, 0d] as double[]
    long[] dh = new long[nbins]
    double dscale = nbins / maxDev
    for (int i = 1; i <= s.getSize(); i++) {
        float[] px = (float[]) s.getPixels(i)
        for (int k = 0; k < px.length; k++) {
            float v = px[k]
            if (v != v) continue
            int b = (int) (Math.abs(v - med) * dscale)
            if (b >= nbins) b = nbins - 1
            dh[b]++
        }
    }
    acc = 0
    int db = 0
    for (int i = 0; i < nbins; i++) {
        acc += dh[i]
        if (acc >= n / 2.0) { db = i; break }
    }
    double mad = maxDev * db / nbins
    return [med, 1.4838 * mad] as double[]
}

// ============================================================================
// 3D FEATURE ENHANCEMENT ALGORITHMS (Native ImgLib2 Refactored)
// ============================================================================

/** Native ImgLib2 3D Difference of Gaussians (DoG) bandpass filter. */
ImagePlus blobDOG(ImagePlus raw, boolean squareRoot, double size) {
    double[] v = voxelSize(raw)
    double s1 = size, s2 = 3 * size
    IJ.log(String.format("  scale %.3f um [%.1f,%.1f,%.1f] px", 
           size, s1 / v[0], s1 / v[1], 2.5 * s1 / v[2]))

    // 1. Wrap raw image directly without memory duplication (zero-copy)
    def imgBase = ImageJFunctions.wrapReal(raw)
    
    // 2. Virtual zero-copy view for square root (if enabled)
    def img = squareRoot ? Converters.convert(imgBase, { src, target -> 
        double val = src.getRealDouble()
        target.setReal(val > 0 ? Math.sqrt(val) : 0.0)
    }, new FloatType()) : imgBase

    // Allocate output buffers for Gaussian scales
    def dog1 = imgBase.factory().create(imgBase)
    def dog2 = imgBase.factory().create(imgBase)

    double[] sigmas1 = [s1 / v[0], s1 / v[1], 2.5 * s1 / v[2]] as double[]
    double[] sigmas2 = [s2 / v[0], s2 / v[1], 2.5 * s2 / v[2]] as double[]

    Gauss3.gauss(sigmas1, Views.extendBorder(img), dog1)
    Gauss3.gauss(sigmas2, Views.extendBorder(img), dog2)

    // Parallel in-place subtraction
    subtractInPlaceImgLib2(dog1, dog2)

    // 3. Converts ImgLib2 VirtualStack -> real, writeable ImageJ1 FloatStack
    ImagePlus out = ImageJFunctions.wrap(dog1, "blobDOG").duplicate()
    out.setCalibration(raw.getCalibration().copy())
    return out
}

/** Native ImgLib2 3D Laplacian of Gaussian (LoG) filter using partial derivatives. */
ImagePlus blobLOG(ImagePlus raw, boolean squareRoot, double size) {
    double[] v = voxelSize(raw)
    
    ImagePlus id1 = duplicateFloat(raw, "id1")
    if (squareRoot) sqrtInPlace(id1)

    // Wrap ImagePlus into ImgLib2 structure
    def img = ImageJFunctions.wrapReal(id1)
    def blurred = img.factory().create(img)
    def laplacianImg = img.factory().create(img)
    def d2 = img.factory().create(img)

    // 1. Parallel Gaussian smoothing
    double[] sigmas = [size / v[0], size / v[1], 2.5 * size / v[2]] as double[]
    Gauss3.gauss(sigmas, Views.extendBorder(img), blurred)

    def extended = Views.extendBorder(blurred)
    int numDimensions = blurred.numDimensions()

    // 2. Sum 2nd partial derivatives across all axes: d^2I / dx^2
    for (int d = 0; d < numDimensions; d++) {
        PartialDerivative.gradientCentralDifference(extended, d2, d, 2)

        LoopBuilder.setImages(laplacianImg, d2).forEachPixel { lap, dev ->
            lap.setReal(lap.getRealDouble() + dev.getRealDouble())
        }
    }

    // 3. Invert sign so bright spots produce positive peaks
    LoopBuilder.setImages(laplacianImg).forEachPixel { pix ->
        pix.setReal(-pix.getRealDouble())
    }

    ImagePlus lap = ImageJFunctions.wrap(laplacianImg, "id1 Laplacian").duplicate()
    dispose(id1)
    lap.setCalibration(raw.getCalibration().copy())
    return lap
}

/** 3D Top-Hat morphology filter (Original - Minimum Rank Filtered). */
ImagePlus blobTH(ImagePlus raw, boolean squareRoot, double size) {
    double[] v = voxelSize(raw)
    double sx = Math.max(3 * size / v[0], 2)
    double sy = Math.max(3 * size / v[1], 2)
    double sz = Math.max(3 * size / v[2], 2)
    ImagePlus id1 = duplicateFloat(raw, "id1")
    if (squareRoot) sqrtInPlace(id1)
    gaussian(id1, 0.75, 0.75, 0.75)
    ImagePlus id2 = rankFilter(id1, false, sx, sy, sz)
    subtractInPlace(id1, id2)
    dispose(id2)
    return id1
}

/** Identifies regional maxima for both 2D images and 3D stacks. */
ImagePlus localMaxima(ImagePlus feat, double size) {
    int slices = feat.getStackSize()
    int w = feat.getWidth()
    int h = feat.getHeight()
    int nPixels = w * h

    ImageStack dst = new ImageStack(w, h)

    if (slices > 1) {
        // 3D Pipeline: 26-connectivity
        ImageStack maxStack = MinimaAndMaxima3D.regionalMaxima(feat.getStack(), 26)
        for (int z = 1; z <= slices; z++) {
            byte[] srcPx = (byte[]) maxStack.getPixels(z)
            float[] floatPx = new float[nPixels]
            for (int i = 0; i < nPixels; i++) {
                floatPx[i] = (srcPx[i] & 0xFF) > 0 ? 1.0f : 0.0f
            }
            dst.addSlice(new FloatProcessor(w, h, floatPx))
        }
    } else {
        // 2D Pipeline: 8-connectivity
        ImageProcessor ip = feat.getProcessor()
        ImageProcessor maxIp = MinimaAndMaxima.regionalMaxima(ip, 8)
        byte[] srcPx = (byte[]) maxIp.getPixels()
        float[] floatPx = new float[nPixels]
        for (int i = 0; i < nPixels; i++) {
            floatPx[i] = (srcPx[i] & 0xFF) > 0 ? 1.0f : 0.0f
        }
        dst.addSlice(new FloatProcessor(w, h, floatPx))
    }

    ImagePlus out = new ImagePlus("localMaxima", dst)
    out.setCalibration(feat.getCalibration().copy())
    return out
}


/**
 * Combined Local Maxima + Thresholding using ImgLib2.
 * Unused currently - need to be compatible with size filtering
 */
ImagePlus combinedLocalMaximaThreshold(ImagePlus feat, double threshold, double spotSize) {
    double[] v = voxelSize(feat)
    int numDims = feat.getStackSize() > 1 ? 3 : 2
    
    // Calculate neighborhood radii in pixels per axis
    int[] span = new int[numDims]
    span[0] = (int) Math.max(1, Math.round(spotSize / v[0]))
    span[1] = (int) Math.max(1, Math.round(spotSize / v[1]))
    if (numDims == 3) {
        span[2] = (int) Math.max(1, Math.round(spotSize / v[2]))
    }

    // Wrap ImagePlus to ImgLib2 view
    def img = ImageJFunctions.wrapReal(feat)
    def outImg = img.factory().create(img)

    // Pad borders with negative infinity so edge peaks are handled correctly
    def extended = Views.extendValue(img, new FloatType(Float.NEGATIVE_INFINITY))

    // Neighborhood shape (skipCenter = true ensures we compare center only against neighbors)
    RectangleShape shape = new RectangleShape(span, true)

    Cursor cursor = Views.interval(img, img).cursor()
    Cursor outCursor = Views.interval(outImg, outImg).cursor()
    def neighborhoodCursor = shape.neighborhoods(extended, img).cursor()

    float thr = (float) threshold

    while (cursor.hasNext()) {
        cursor.fwd()
        outCursor.fwd()
        neighborhoodCursor.fwd()

        float val = cursor.get().getRealFloat()

        // 1. FAST SHORT-CIRCUIT: Skip neighbor check entirely if below threshold
        if (val < thr) {
            outCursor.get().setReal(0.0f)
            continue
        }

        // 2. NEIGHBORHOOD CHECK: Only runs on candidate foreground voxels
        boolean isMax = true
        Neighborhood nh = neighborhoodCursor.get()
        for (FloatType neighbor : nh) {
            if (neighbor.getRealFloat() >= val) { // Strict local peak condition
                isMax = false
                break
            }
        }

        outCursor.get().setReal(isMax ? 1.0f : 0.0f)
    }

    // Convert ImgLib2 result stack back to ImageJ FloatStack
    ImagePlus out = ImageJFunctions.wrap(outImg, "spots").duplicate()
    out.setCalibration(feat.getCalibration().copy())
    return out
}

/** Applies adaptive MAD-based or fixed cutoff thresholding to feature map. */
void thresholdAuto(ImagePlus imp, double pfa, boolean adaptiveThreshold) {
    double[] mm = minMax(imp)
    double threshold
    if (adaptiveThreshold) {
        double[] stats = robustStats(imp, mm[0], mm[1])
        double lambda = -2 * normppf(Math.pow(10, -pfa))
        threshold = stats[0] + lambda * stats[1]
        IJ.log("    adaptive threshold")
        IJ.log("    min: " + mm[0] + ", max: " + mm[1])
        IJ.log("    mean " + stats[0] + ", std: " + stats[1])
        IJ.log("    specificity: " + pfa + ", quantile: " + lambda + ", threshold: " + threshold)
    } else {
        IJ.log("    fixed threshold")
        IJ.log("    min: " + mm[0] + ", max: " + mm[1])
        IJ.log("    threshold: " + pfa)
        threshold = pfa
    }
    binarizeGE(imp, threshold)
}

// ============================================================================
// MORPHOLIBJ SEGMENTATION & VOLUME FILTERING
// ============================================================================

/** Inverts image stack gray values in place. */
void invertInPlace(ImagePlus imp) {
    ImageStack s = imp.getStack()
    for (int i = 1; i <= s.getSize(); i++) s.getProcessor(i).invert()
}

/** Converts mask volume to 32-bit FloatProcessor binary stack [0.0, 1.0]. */
ImagePlus toBinaryFloat(ImagePlus imp, String title) {
    ImageStack s = imp.getStack()
    ImageStack d = new ImageStack(imp.getWidth(), imp.getHeight())
    for (int i = 1; i <= s.getSize(); i++) {
        ImageProcessor ip = s.getProcessor(i)
        int n = imp.getWidth() * imp.getHeight()
        float[] px = new float[n]
        for (int k = 0; k < n; k++) px[k] = ip.getf(k) > 0 ? 1f : 0f
        d.addSlice(new FloatProcessor(imp.getWidth(), imp.getHeight(), px))
    }
    ImagePlus out = new ImagePlus(title, d)
    out.setCalibration(imp.getCalibration().copy())
    return out
}

/** Segments spot boundaries via 3D MorphoLibJ Watershed and removes objects > maxSize. */
ImagePlus filterOutMasks(ImagePlus input, ImagePlus marker, double maxSize) {
    multiplyInPlace(marker, input)

    double[] v = voxelSize(input)
    double maxVol = 4.0 / 3.0 * Math.PI * (maxSize / v[0]) * (maxSize / v[1]) * (maxSize / v[2])
    IJ.log("    size filtering: " + maxSize + " um radius / " + maxVol + " voxels volume")

    ImagePlus mask8 = BinaryImages.binarize(input)
    mask8.setCalibration(input.getCalibration().copy())
    ImagePlus marker8 = BinaryImages.binarize(marker)
    marker8.setCalibration(marker.getCalibration().copy())
    ImagePlus markerLbl = BinaryImages.componentsLabeling(marker8, 6, 16)
    IJ.log("    labeled " + LabelImages.findAllLabels(markerLbl).length + " markers")

    ImagePlus dist = BinaryImages.distanceMap(mask8)
    invertInPlace(dist)

    ImagePlus water = Watershed.computeWatershed(dist, markerLbl, mask8, 6, false)
    int[] labelsBefore = LabelImages.findAllLabels(water)
    IJ.log("    number of label in the image before filtering: " + labelsBefore.length)

    int[] voxelCounts = LabelImages.voxelCount(water.getStack(), labelsBefore)
    ImageStack filteredStack = water.getStack().duplicate()

    for (int z = 1; z <= filteredStack.getSize(); z++) {
        ImageProcessor ip = filteredStack.getProcessor(z)
        for (int p = 0; p < ip.getPixelCount(); p++) {
            int label = ip.get(p)
            boolean keep = false
            for (int i = 0; i < labelsBefore.length; i++) {
                if (label == labelsBefore[i] && voxelCounts[i] <= maxVol) {
                    keep = true
                    break
                }
            }
            if (!keep) ip.set(p, 0)
        }
    }
    ImagePlus output = new ImagePlus("filtered", filteredStack)
    output.setCalibration(water.getCalibration().copy())
    return output
}

/** Localizes all voxels above threshold back into a coordinate list. */
List localizeSpots(ImagePlus imp, int channelIdx) {
    List coords = []
    int w = imp.getWidth(), h = imp.getHeight()
    ImageStack s = imp.getStack()
    for (int z = 1; z <= s.getSize(); z++) {
        float[] px = (float[]) s.getPixels(z)
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                if (px[y * w + x] >= 0.5f) coords << ([channelIdx, x, y, z] as double[])
            }
        }
    }
    return coords
}

/**
 * Samples intensity along a line centered at (x, y, z) with directional vector (u, v, w).
 * @return 2D double array [4][n] containing [x_coords, y_coords, z_coords, intensities]
 */
def lineIntensity(ImagePlus imp, double x, double y, double z, double u, double v, double w, int n) {
    int width = imp.getWidth(), height = imp.getHeight(), slices = imp.getStackSize()
    def img = ImageJFunctions.wrapReal(imp)
    def ra = Views.interpolate(Views.extendBorder(img), new NLinearInterpolatorFactory()).realRandomAccess()

    def a = new double[4][n]
    for (int i = 0; i < n; i++) {
        double t = (double) i / (n - 1) - 0.5
        double xi = Math.min(width - 1, Math.max(x + t * u, 0))
        double yi = Math.min(height - 1, Math.max(y + t * v, 0))
        double zi = Math.min(slices - 1, Math.max(z - 1 + t * w, 0))

        ra.setPosition([xi, yi, zi] as double[])

        a[0][i] = xi
        a[1][i] = yi
        a[2][i] = zi + 1
        a[3][i] = ra.get().getRealDouble()
    }
    return a
}

/** Fits a 1D Gaussian curve to sampled intensity profile and returns center parameter. */
double fitCenter(double[] xs, double[] ys) {
    CurveFitter cf = new CurveFitter(xs, ys)
    cf.doFit(CurveFitter.GAUSSIAN)
    return cf.getParams()[2]
}

/** Refines 3D position using 1D orthogonal Gaussian profile fitting along X, Y, Z. */
double[] gaussianRefine(ImagePlus imp, double[] p) {
    int n = 5
    int slices = imp.getStackSize()
    double x = p[0], y = p[1], z = p[2]
    try {
        for (int iter = 0; iter < 5; iter++) {
            def a = lineIntensity(imp, x, y, z, n, 0, 0, n)
            x = fitCenter(a[0], a[3])
            a = lineIntensity(imp, x, y, z, 0, n, 0, n)
            y = fitCenter(a[1], a[3])
            if (slices > n) {
                a = lineIntensity(imp, x, y, z, 0, 0, n, n)
                z = fitCenter(a[2], a[3])
            }
            if (Double.isNaN(x) || Double.isNaN(y) || Double.isNaN(z)) return p
        }
    } catch (Exception e) {
        return p
    }
    if (Math.abs(x - p[0]) < 1 && Math.abs(y - p[1]) < 1 && Math.abs(z - p[2]) < 1) {
        return [x, y, z] as double[]
    }
    return p
}

/** Subpixel refinement. */
void refineLocalization(ImagePlus raw, List coords) {
    coords.each { double[] c ->
        double[] p = [c[1], c[2], c[3]] as double[]
        p = gaussianRefine(raw, p)
        c[1] = p[0] 
        c[2] = p[1] 
        c[3] = p[2] 
    }
}

/** Converts voxel coordinates to calibrated micron space with optional subpixel refinement. */
void scaleCoordinatesToMicrons(ImagePlus raw, List coords) {	
	double[] d = voxelSize(raw)
	coords.each { double[] c ->
        double[] p = [c[1], c[2], c[3]] as double[]
        c[1] = p[0] * d[0]
        c[2] = p[1] * d[1]
        c[3] = p[2] * d[2]
    }
}

/** Single-channel spot detection workflow execution block. */
List detect3DSpots(ImagePlus imp, List channels, String feat, List sizes,
        List pfas, int channelIdx, boolean subpix, ImagePlus maskImp, boolean adaptiveThr,
        List maxSizes) {
    int channel = channels[channelIdx]
    double size = sizes[channelIdx]
    double pfa = pfas[channelIdx]
    double maxSize = maxSizes[channelIdx]

    IJ.log(" - detect spots in channel " + channel)
    
    ImagePlus raw = extractChannelFloat(imp, channel)
    ImagePlus id1
    // Step 1 - Compute a blob feature
    switch (feat.toUpperCase()) {
        case "DOG": id1 = blobDOG(raw, true, size); break
        case "LOG": id1 = blobLOG(raw, true, size); break
        case "TOP HAT": id1 = blobTH(raw, true, size); break
        default: id1 = duplicateFloat(raw, "id1")
    }
    
	// Step 2 - Non-maximuum suppression
    ImagePlus id3 = localMaxima(id1, size)
    
    // Step 3 - Apply predefined masking on the NMS
    if (maskImp != null) multiplyInPlace(id3, maskImp)
    
	// Step 4 - Threshold the score in place
    thresholdAuto(id1, pfa, adaptiveThr)
    
	// Step 5 - Filter out the segmented region above maxSize
    if (maxSize > 0) {
        ImagePlus filtered = filterOutMasks(id1, id3, maxSize)
        dispose(id1)
        id1 = filtered
    }
    
	//  Step 6 - Multiply/And operator of thresholded map id1 and NMS id3
    multiplyInPlace(id1, id3)
    
	//  Step 7 - Collect pixel locations
    IJ.log("   localize spots")
    List coords = localizeSpots(id1, channelIdx)

	if (subpix) {
	    IJ.log("   refine spots")
	    refineLocalization(raw, coords)
	}
	scaleCoordinatesToMicrons(raw, coords)
	

    [id1, id3, raw].each { dispose(it) }
    IJ.log("   number of spots: " + coords.size())
    return coords
}

/** Iterates spot detection across all user-selected channels. */
List detectSpotsInAllChannels(ImagePlus imp, List channels, String feat,
        List specificity, List spotSizes, ImagePlus maskImp, boolean adaptiveThr,
        List maxSizes, boolean subpix) {
    IJ.log("Detect 3D spots in image " + imp.getTitle())
    List coords = []
    try {
        for (int idx = 0; idx < channels.size(); idx++) {
            coords.addAll(detect3DSpots(imp, channels, feat, spotSizes, specificity, idx, subpix,
                maskImp, adaptiveThr, maxSizes))
        }
    } finally {
        print("success")
    }
    return coords
}

/** Prepares binary ROI mask volume from a specified channel. */
ImagePlus getMask(ImagePlus imp, int channel) {
    IJ.log("Using channel " + channel + " as segmentation mask.")
    ImagePlus m = extractChannelFloat(imp, channel)
    binarizeGE(m, 1.0)
    m.setTitle("segmentation")
    return m
}

// ============================================================================
// DATA SERIALIZATION & IJ RESULTS TABLE UTILITIES
// ============================================================================

/** Converts list of spot coordinates into ImageJ ResultsTable. */
ResultsTable coords2Table(List coords, List channels, List sizes) {
    ResultsTable rt = new ResultsTable()
    coords.each { double[] c ->
        int ci = (int) c[0]
        rt.incrementCounter()
        rt.addValue("Channel", channels[ci])
        rt.addValue("X [um]", c[1])
        rt.addValue("Y [um]", c[2])
        rt.addValue("Z [um]", c[3])
        rt.addValue("Size [um]", sizes[ci])
    }
    return rt
}

/** Reads spot coordinates from an external CSV file. */
List loadCoordsTable(String path, List channels) {
    ResultsTable rt = ResultsTable.open(path)
    List coords = []
    for (int i = 0; i < rt.size(); i++) {
        int idx = channels.indexOf((int) rt.getValue("Channel", i))
        if (idx < 0) continue
        coords << ([idx, rt.getValue("X [um]", i), rt.getValue("Y [um]", i), rt.getValue("Z [um]", i)] as double[])
    }
    return coords
}

// ============================================================================
// SPATIAL COLOCALIZATION ENGINE (ImgLib2 KD-Tree)
// ============================================================================

/**
 * Builds 3D KD-Tree across all spot coordinates and calculates neighborhood channel occupancy codes.
 *
 * @param coords List of spot detections [channelIdx, x_um, y_um, z_um]
 * @param nc Total number of analyzed channels
 * @param dmaxUm Spatial proximity cutoff radius in micrometers
 * @return Flat integer array representing channel occupancy counts per spot
 */
int[] computeCodes(List coords, int nc, double dmaxUm) {
    int N = coords.size()
    if (N == 0) return new int[0]

    List points = new ArrayList<>(N)
    List dummyValues = new ArrayList<>(N)

    for (int i = 0; i < N; i++) {
        double[] c = coords[i]
        points.add(new RealPoint(c[1], c[2], c[3]))
        dummyValues.add(new DoubleType(c[0]))
    }

    KDTree tree = new KDTree<>(dummyValues, points)
    RadiusNeighborSearchOnKDTree search = new RadiusNeighborSearchOnKDTree<>(tree)

    int[] codes = new int[nc * N]
    for (int i = 0; i < N; i++) {
        RealPoint center = points.get(i)
        search.search(center, dmaxUm, true)

        for (int k = 0; k < search.numNeighbors(); k++) {
            int neighborChannel = (int) search.getSampler(k).get().getRealDouble()
            codes[neighborChannel + i * nc]++
        }
    }
    return codes
}

/** Generates binary matrix representing power set (2^n combinations) of channel sets. */
int[] powerSet(int n) {
    int m = 1 << n
    int[] set = new int[n * m]
    for (int i = 0; i < m; i++) {
        String bin = Integer.toBinaryString(i).padLeft(n, "0")
        for (int j = 0; j < n; j++) set[j + n * i] = bin.charAt(j) == ('1' as char) ? 1 : 0
    }
    return set
}

/** Evaluates occupancy counts against power-set channel combinations. */
double[] countsCodeBySet(double[] codes, int[] sets, int n) {
    int m = 1 << n
    int N = codes.length / n
    double[] counts = new double[N * m]
    for (int i = 0; i < N; i++) {
        for (int k = 1; k < m; k++) {
            double acc1 = 0, acc2 = 0
            boolean test = true
            for (int j = 0; j < n; j++) {
                acc1 += codes[j + n * i] * sets[j + n * k]
                acc2 += sets[j + n * k]
                if (sets[j + n * k] == 1) test = test && (codes[j + n * i] > 0)
            }
            if (test) counts[i + N * k] = acc1 / acc2
        }
    }
    return counts
}

/** Aggregates colocalized spot counts across channel power set. */
double[] aggregateCountsPerSet(double[] counts, int n) {
    int m = 1 << n
    int N = counts.length / m
    double[] agg = new double[m]
    for (int k = 0; k < m; k++) {
        for (int i = 0; i < N; i++) if (counts[i + N * k] > 0) agg[k]++
    }
    return agg
}

/** Checks bitwise subset inclusion between two channel binary sets. */
boolean testIntersection(int[] code1, int[] code2) {
    if (code2.length == 0) return false
    for (int i = 0; i < code1.length; i++) {
        if (code2[i] == 1 && code1[i] != 1) return false
    }
    return true
}

/** Applies inclusion-exclusion correction to remove higher-order colocalization double-counting. */
double[] correctAggCounts(double[] agg, int[] sets, int n) {
    double[] agg2 = (double[]) agg.clone()
    int m = 1 << n
    for (int j = m - 1; j > 1; j--) {
        int[] a = Arrays.copyOfRange(sets, j * n, (j + 1) * n)
        for (int i = j - 1; i > 0; i--) {
            int[] b = Arrays.copyOfRange(sets, i * n, (i + 1) * n)
            if (testIntersection(a, b)) agg2[i] = agg2[i] - agg2[j]
        }
    }
    return agg2
}

/** Normalizes aggregate colocalization counts by channel set cardinality. */
double[] divideBySetCardinality(double[] agg, int[] sets, int n) {
    double[] agg2 = (double[]) agg.clone()
    for (int k = 1; k < agg.length; k++) {
        int card = 0
        for (int j = 0; j < n; j++) card += sets[j + n * k]
        agg2[k] = agg[k] / card
    }
    agg2[0] = 0
    return agg2
}

// ============================================================================
// UI & SUMMARY REPORTING UTILITIES
// ============================================================================

/** Constructs channel combination title string (e.g., "Ch1Ch2"). */
String getSetName(int[] set, List channels) {
    String str = ""
    for (int j = 0; j < channels.size(); j++) if (set[j] > 0) str += "Ch" + channels[j]
    return str
}

/** Sorts channel set indices by set cardinality for clean tabular presentation. */
List rankSetByCardinality(int[] sets, int n) {
    int m = 1 << n
    List order = new ArrayList(m)
    for (int k = 0; k < m; k++) order.add(k)
    order.sort { a, b ->
        int cardA = 0
        int cardB = 0
        for (int j = 0; j < n; j++) {
            cardA += sets[j + n * a]
            cardB += sets[j + n * b]
        }
        Integer.compare(cardA, cardB)
    }
    return order
}

/** Renders per-spot channel code table in ImageJ GUI. */
void showCodes(String tbl, int[] codes, List channels) {
    int n = channels.size()
    int N = codes.length / n
    ResultsTable rt = new ResultsTable()
    for (int i = 0; i < N; i++) rt.incrementCounter()
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < n; j++) rt.setValue(String.valueOf(channels[j]), i, codes[j + n * i])
    }
    rt.show(tbl)
}

/** Renders per-spot combination membership counts in ImageJ GUI. */
void showCounts(String tbl, double[] counts, int[] sets, List channels) {
    int n = channels.size()
    int m = 1 << n
    int N = counts.length / m
    List rank = rankSetByCardinality(sets, n)
    ResultsTable rt = new ResultsTable()
    for (int i = 0; i < N; i++) rt.incrementCounter()
    for (int k = 1; k < m; k++) {
        int ko = rank[k]
        String col = getSetName(Arrays.copyOfRange(sets, n * ko, n * (ko + 1)), channels)
        for (int i = 0; i < N; i++) rt.setValue(col, i, counts[i + N * ko])
    }
    rt.show(tbl)
}

/** Appends final colocalization statistics row to persistent summary table 'MultiChannelsColocalisation.csv'. */
void appendRecord(double[] agg, int[] sets, List channels, String feat,
        List specificity, List spotSize, String name) {
    String tbl = "MultiChannelsColocalisation.csv"
    Frame win = WindowManager.getFrame(tbl)
    ResultsTable rt = (win instanceof TextWindow) ? ((TextWindow) win).getTextPanel().getResultsTable() : null
    if (rt == null) rt = new ResultsTable()

    int n = channels.size()
    int m = 1 << n
    rt.incrementCounter()
    rt.addValue("Name", name)
    List rank = rankSetByCardinality(sets, n)
    double total = 0
    for (int k = 1; k < m; k++) {
        int ko = rank[k]
        String col = getSetName(Arrays.copyOfRange(sets, n * ko, n * (ko + 1)), channels)
        double val = (double) Math.round(agg[ko])
        rt.addValue(col, val)
        total += val
    }
    rt.addValue("Number", total)
    rt.addValue("Feature", feat)
    for (int i = 0; i < n; i++) {
        rt.addValue("size Ch" + channels[i], spotSize[i])
        rt.addValue("specificity Ch" + channels[i], specificity[i])
    }
    rt.show(tbl)
}

// ============================================================================
// SYNTHETIC TEST BENCH GENERATOR
// ============================================================================

/** Synthesizes synthetic 3D hyperstack with ground truth spot distributions and Gaussian noise. */
ImagePlus generateTestImage(List channels, String recordName) {
    closeAllImages()
    int w = 200, h = 200, d = 100, n = 100, b = 10
    int c = channels.size()
    IJ.log("Generating " + n + " candidate locations...")
    Random rnd = new Random()

    ImageStack[] chs = new ImageStack[c]
    for (int j = 0; j < c; j++) {
        chs[j] = new ImageStack(w, h)
        for (int z = 0; z < d; z++) {
            ShortProcessor sp = new ShortProcessor(w, h)
            Arrays.fill((short[]) sp.getPixels(), (short) 100)
            chs[j].addSlice(sp)
        }
    }
	IJ.log("Generating image...")
    int m = 1 << c
    double[] agg0 = new double[m]
    List codesList = []
    for (int i = 0; i < n; i++) {
        int x = b + (int) Math.round((w - 2 * b) * rnd.nextDouble())
        int y = b + (int) Math.round((h - 2 * b) * rnd.nextDouble())
        int z = b + (int) Math.round((d - 2 * b) * rnd.nextDouble())
        int[] values = new int[c]
        for (int j = 0; j < c; j++) {
            if (rnd.nextDouble() > 0.2) {
                ((short[]) chs[j].getPixels(z))[y * w + x] = (short) 10000
                values[j] = 1
            }
        }
        for (int j0 = 0; j0 < c; j0++) if (values[j0] == 1) codesList.addAll(values.toList())
        agg0[Integer.parseInt(values.toList().join(""), 2)]++
    }
    
	IJ.log("Generating blur and noise...")
    double sigma = 0.2, dxy = 0.11, dz = 0.33
    for (int j = 0; j < c; j++) {
        ImageStack st = Filters3D.filter(chs[j], Filters3D.MAX, 1f, 1f, 1f)
        ImagePlus tmp = new ImagePlus("ch" + (j + 1), st)
        gaussian(tmp, sigma / dxy, sigma / dxy, 3 * sigma / dz)
        st = tmp.getStack()
        for (int z = 1; z <= st.getSize(); z++) st.getProcessor(z).noise(25.0)
        chs[j] = st
    }
	
	IJ.log("Assembling stack...")
    ImageStack all = new ImageStack(w, h)
    for (int z = 1; z <= d; z++) for (int j = 0; j < c; j++) all.addSlice(chs[j].getProcessor(z))
    ImagePlus imp = new ImagePlus("Test Image", all)
    imp.setDimensions(c, d, 1)
    if (c > 1) imp = new CompositeImage(imp, IJ.GRAYSCALE)
    imp.setOpenAsHyperStack(true)
    
    Calibration cal = imp.getCalibration()
    cal.pixelWidth = dxy; cal.pixelHeight = dxy; cal.pixelDepth = dz; cal.setUnit("um")
    imp.setSlice((int) Math.round(d / 2.0))
    for (int i = 1; i <= c; i++) {
        imp.setC(i)
        imp.setDisplayRange(50, 1000)
    }
    
	IJ.log("Testing aggregation...")
    double[] codes = codesList.collect { it as double } as double[]
    int[] sets = powerSet(c)
    double[] counts = countsCodeBySet(codes, sets, c)
    double[] agg = aggregateCountsPerSet(counts, c)
    double[] agg4 = divideBySetCardinality(correctAggCounts(agg, sets, c), sets, c)
    
    IJ.log("Storing results..")
    List zeros = [0d] * c

    appendRecord(agg4, sets, channels, "ground truth", zeros, zeros, recordName)
    appendRecord(agg0, sets, channels, "ground truth*", zeros, zeros, recordName)
    return imp
}

// ============================================================================
// VISUALIZATION & OVERLAY EXPORT
// ============================================================================

/** Derives channel overlay color from display LUT or defaults to yellow. */
Color colorFromLut(ImagePlus imp, int channel) {
    try {
        LUT[] luts = imp.getLuts()
        if (luts != null && channel >= 1 && channel <= luts.length) {
            LUT l = luts[channel - 1]
            return new Color((int) (0.2 * 255 + 0.8 * l.getRed(255)),
                (int) (0.2 * 255 + 0.8 * l.getGreen(255)),
                (int) (0.2 * 255 + 0.8 * l.getBlue(255)))
        }
    } catch (Exception e) { /* fall through */ }
    return Color.YELLOW
}

/** Generates 3D circular ROI overlay markers on localized spot coordinates. */
void addOverlay(ImagePlus imp, List coords, List channels, List spotSizes) {
    double[] d = voxelSize(imp)
    int nz = imp.getNSlices()
    Overlay ov = new Overlay()
    coords.each { double[] c ->
        int ci = (int) c[0]
        double x = c[1] / d[0], y = c[2] / d[1]
        int z = Math.max(1, Math.min(nz, (int) Math.round(c[3] / d[2])))
        double R = 2 * spotSizes[ci] / d[0]
        OvalRoi roi = new OvalRoi(x - R, y - R, 2 * R + 1, 2 * R + 1)
        roi.setStrokeColor(colorFromLut(imp, channels[ci]))
        roi.setPosition(channels[ci], nz > 1 ? z : 1, 0)
        ov.add(roi)
    }
    imp.setOverlay(ov)
}

/** Exports detected spot coordinate point table to external CSV file. */
void exportCoordsToCSV(String mode, List coords, List channels, List spotSizes,
        ImagePlus imp, String recordName, boolean save) {
    if (mode == "csv" || closeOnExit) return
    ResultsTable rt = coords2Table(coords, channels, spotSizes)
    if (save) {
        File f = new File(recordName)
        File dir = f.isAbsolute() ? f.getParentFile() : null
        if (dir == null && imp != null && imp.getOriginalFileInfo() != null && imp.getOriginalFileInfo().directory) {
            dir = new File(imp.getOriginalFileInfo().directory)
        }
        if (dir == null) dir = new File(System.getProperty("user.home"))
        File out = new File(dir, stripExtension(f.getName()) + "-points.csv")
        rt.saveAs(out.getAbsolutePath())
        IJ.log("Coordinates saved to " + out.getAbsolutePath())
    }
    rt.show("points.csv")
}

/** Generates Maximum Intensity Projection (MIP) and displays ROI overlay. */
void zProjectAndShowROIs(ImagePlus imp, List coords, List channels, List spotSizes) {
    imp.deleteRoi()
    ImagePlus proj = imp
    if (imp.getNSlices() > 1) {
        proj = ZProjector.run(imp, "max all")
        proj.show()
    }
    addOverlay(proj, coords, channels, spotSizes)
    if (proj.isComposite()) proj.setDisplayMode(IJ.COMPOSITE)
    for (int c = 1; c <= proj.getNChannels(); c++) {
        proj.setC(c)
        IJ.run(proj, "Enhance Contrast", "saturated=0.35")
    }
}

/** Load the image using Bio-Formats **/
ImagePlus openImage(String path) {
    ImagePlus[] imps = BF.openImagePlus(path)
    if (imps != null && imps.length > 0) {
        return imps[0]
    } else {
        IJ.error("Error", "Failed to open image via Bio-Formats: " + path)
        return null
    }
}

// ============================================================================
// STATISTICAL ANALYSIS
// ============================================================================

/**
 * Calculates Poisson p-values for all power-set channel combinations against CSR
 * Very crude approximation
 * 
 */
//Map computeCSRPvalues(
//	double[] channelCounts, 
//    double[] observedAggCounts, 
//    int[] powerSets, 
//    int nc, 
//    double dmax, 
//    double totalVolumeUm, 
//    boolean is3D) {
//    Map pValues = [:]
//    if (totalVolumeUm <= 0) return pValues
//
//    // Search ball volume
//    double vR = is3D ? (4.0 / 3.0) * Math.PI * Math.pow(dmax, 3) 
//                     : Math.PI * Math.pow(dmax, 2)
//                     
//    for (int k = 0; k < powerSets.size(); k++) {
//    	
//    }
//                     
//    retrun pValues
//        	
//}

/**
 * Calculates Poisson p-values for all power-set channel combinations.
 * 
 * @param channelCounts Array containing total spot counts per channel [N_ch1, N_ch2, ...]
 * @param observedAggCounts Array of observed aggregate spot counts per set (from aggregateCountsPerSet)
 * @param powerSets Binary power-set matrix (from powerSet(nc))
 * @param nc Number of channels
 * @param dmax Proximity radius in microns
 * @param totalVolumeUm Physical volume of mask or image bounding box (um^3 or um^2)
 * @param is3D True if 3D stack, false if 2D slice
 * @return Map of set index to p-value
 */
Map computeMultiChannelPoissonPValues(
        double[] channelCounts, 
        double[] observedAggCounts, 
        int[] powerSets, 
        int nc, 
        double dmax, 
        double totalVolumeUm, 
        boolean is3D) {

    Map pValues = [:]
    if (totalVolumeUm <= 0) return pValues

    // Search ball volume
    double vR = is3D ? (4.0 / 3.0) * Math.PI * Math.pow(dmax, 3) 
                     : Math.PI * Math.pow(dmax, 2)

    int numSets = 1 << nc

    // Iterate across power-set combinations (skip empty set k=0 and singletons card=1)
    for (int k = 1; k < numSets; k++) {
        List memberChannels = []
        List nonMemberChannels = []

        for (int j = 0; j < nc; j++) {
            if (powerSets[j + nc * k] == 1) {
                memberChannels.add(j)
            } else {
                nonMemberChannels.add(j)
            }
        }

        // Only process multi-channel combinations (2 or more channels)
        if (memberChannels.size() < 2) {
            pValues[k] = 1.0
            continue
        }

        // Select reference channel with minimum spot count
        int refChan = memberChannels.min { channelCounts[it] }
        int nRef = channelCounts[refChan]

        if (nRef == 0) {
            pValues[k] = 1.0
            continue
        }

        // Calculate joint probability p_joint (Inclusive model)
        double pJoint = 1.0
        for (int targetChan : memberChannels) {
            if (targetChan == refChan) continue
            double lambdaJ = (double) channelCounts[targetChan] / totalVolumeUm
            double pJ = 1.0 - Math.exp(-lambdaJ * vR)
            pJoint *= pJ
        }

        double lambdaTotal = nRef * pJoint
        int kObs = (int) Math.round(observedAggCounts[k])

        // Poisson CDF survival function: P(K >= kObs) = 1 - CDF(kObs - 1)
        double sum = 0.0
        double term = Math.exp(-lambdaTotal)
        sum += term

        for (int i = 1; i < kObs; i++) {
            term *= (lambdaTotal / i)
            sum += term
        }

        double pVal = 1.0 - sum
        pValues[k] = Math.max(0.0, Math.min(1.0, pVal))
    }
	print(pValues)
    return pValues
}

// ============================================================================
// MAIN WORKFLOW 
// ============================================================================

/** Main entry point for workflow execution. */
void runColocalization() {
    IJ.log("Starting Colocalization Workflow...")
    long t0 = System.currentTimeMillis()

    inp = inputPath.getAbsolutePath()
    mode = getMode(inp)
    IJ.log("Execution Mode: " + mode)

    List channels = parseCSVInt(channelsStr)
    List specificity = parseCSVFloat(specificityStr)
    List spotSizes = parseCSVFloat(spotSizeStr)
    List maxSizes = parseCSVFloat(maxSpotSizeStr)
    Integer maskChannel = parseMaskStr(maskStr)
    boolean adaptiveThr = adaptive
    String recordName = inp

    if (mode == "test") {
        recordName = "Test Image.tif"
        channels = [1, 2, 3]
        spotSizes = [0.2d] * 3
        specificity = [3d] * 3
        adaptiveThr = true
        maskChannel = null
    }

    int nc = channels.size()
    if ([specificity, spotSizes, maxSizes].any { it.size() < nc }) {
        IJ.error("Multichannel spot colocalization",
            "Spot size, max spot size and specificity need one value per channel (" + nc + ").")
        return
    }

    ImagePlus imp = null
    List coords
    if (mode == "test") {
        IJ.log("Generating synthetic test dataset...")
        imp = generateTestImage(channels, recordName)
        IJ.log("Detecting...")
        coords = detectSpotsInAllChannels(imp, channels, feature, specificity, spotSizes, null,
            adaptiveThr, maxSizes, subpixel)
    } else if (mode == "csv") {
        IJ.log("Loading csv " + inp + " table")
        coords = loadCoordsTable(inp, channels)
    } else {
        if (mode == "image") {
            if (WindowManager.getImageCount() == 0) {
                String p = IJ.getFilePath("Please select a file")
                if (p == null) return
                imp = IJ.openImage(p)
            } else {
                imp = WindowManager.getCurrentImage()
            }
            IJ.log("[Using active image " + imp.getTitle() + "]")
            recordName = imp.getTitle()
        } else {
            IJ.log("Opening image " + inp)
            imp = openImage(inp)
            if (imp == null) return
        }
        ImagePlus maskImp = maskChannel != null ? getMask(imp, maskChannel) : null
        IJ.log("Detecting...")
        coords = detectSpotsInAllChannels(imp, channels, feature, specificity, spotSizes, maskImp,
            adaptiveThr, maxSizes, subpixel)
        dispose(maskImp)
    }

    IJ.log("Counting and aggregating...")
    int[] icodes = computeCodes(coords, nc, dmax)
    showCodes("codes.csv", icodes, channels)
    double[] codes = icodes.collect { it as double } as double[]
    int[] sets = powerSet(nc)
    double[] counts = countsCodeBySet(codes, sets, nc)
    showCounts("counts.csv", counts, sets, channels)
    double[] agg = aggregateCountsPerSet(counts, nc)
    double[] agg2 = divideBySetCardinality(correctAggCounts(agg, sets, nc), sets, nc)
    
	IJ.log("Calculating csr pvalues...")
	double volumeUm = computePhysicalVolume(imp)
	boolean is3D = imp.nSlices > 1 ? true : false
    pvalues = computeMultiChannelPoissonPValues(agg, agg2, sets, nc, dmax, volumeUm, is3D)
    print(pvalues)
    
    appendRecord(agg2, sets, channels, feature, specificity, spotSizes, recordName)
    
    if (imp != null) imp.show()
    
    if (imp != null && !closeOnExit) addOverlay(imp, coords, channels, spotSizes)
    exportCoordsToCSV(mode, coords, channels, spotSizes, imp, recordName, saveCoordinates)
    if (doZProject && imp != null && !closeOnExit) zProjectAndShowROIs(imp, coords, channels, spotSizes)

    long t1 = System.currentTimeMillis()
    System.gc()
    IJ.log("Finished in " + (t1 - t0) / 1000.0 + " seconds.")
    if (closeOnExit) closeAllImages()
}

// Execute pipeline
runColocalization()