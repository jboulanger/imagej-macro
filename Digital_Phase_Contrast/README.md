# Digital Phase Contrast

Convert a brightfield image pair or stack to a quantitative phase image.

The macro can take as input 3D stack or 3D timelapses hyperstacks.

![image](https://user-images.githubusercontent.com/3415561/219726074-31f33b51-7a36-4b96-a0ca-0618d68d31a3.png)

The two extreme images of the brightfield stack are subtracted from one another and a Poisson solver is used on the difference image to retreive the phase. A post processing finally remove background and normalize images to 16 bits.


