# Beads_FWHM

This macro automates the esimtation the full width half maximum (FWHM) of imaged fluorescent beads.

Beads are detected automatically and several line profile extracted. For each line profile, a 1D gaussian function is fitted. From the average standard deviation of the Gaussian fit several resolution criterion are deducted. The relationship was obtined by fitting an Airy function with a Gaussian function (see airy_gaussian_approx.m script). The appearance of the beads will also depends on their physical size. We found that approximating the top hat profile of width b of a bead and a Gaussian function of standard deviation b/2 was allowing to compensate the bead size for small beads.

## Tutorial
- Open an image. If no image are opened a test image with 25 regularly spaced beads is generated.
![image](https://user-images.githubusercontent.com/3415561/122259643-4c139700-ceca-11eb-8d2c-795bba9b4560.png)
- Run the macro and fill in the parameters
![image](https://user-images.githubusercontent.com/3415561/122259531-31412280-ceca-11eb-83e2-aa92b5a7086a.png)
- The macro generate a table with the parameter of each object:
![image](https://user-images.githubusercontent.com/3415561/122259803-75ccbe00-ceca-11eb-86e0-88e492d81a6b.png)




