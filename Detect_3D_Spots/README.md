# Detect_3D_Spots

Detect spots in 3D hyperstacks using a Difference of Gaussian filter. Count the number of spots within user defined regions of interest. Optionally, spots are rendered on the image.

If no images are opened, a test image is generated.

![image](https://user-images.githubusercontent.com/3415561/117969900-ea3a9d00-b31f-11eb-9d28-65b599a813f0.png)


## Tutorial 
1. Open a 3D hyperstack or close all images to generate a test image
2. Add some region of interest to the ROI manager (press t)
3. Open the macro and press run.
4. Set the parameters:
  - Channel: set the channel that will be process for spot detection
  - Spot size: set the approximate size in pixel of the spot to detect 
  - Probability of false alarm: set the probability of false detection (log10 scale)
  - Render: Choose if you want to render a control image
  - Rendering spot size: set the size in pixels of the rendered spots.

![image](https://user-images.githubusercontent.com/3415561/117969971-ffafc700-b31f-11eb-9b9f-579b4e784820.png)

5.Inspect the results: table and rendered image.
