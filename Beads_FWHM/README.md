# Beads_FWHM

This macro estimate the full width half maximum (FWHM) of imaged fluorescent beads by fitting a Gaussian function to the bead images.

From the standard deviation of the Gaussian fit several resolution criterion are deducted. The relationship was obtined by fitting an Airy function with a Gaussian function. The appearance of the beads will also depends on their physical size. We found that approximating the top hat profile of width b of a bead and a Gaussian function of standard deviation b/2 was allowing to compensate the bead size for small beads.

If no image are opened a test image with 25 regularly spaced beads is generated.

## Tutorial

- Open an image
- Run the macro and fill in the parameters

