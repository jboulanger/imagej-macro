# Volumic intensity analysis

This macro computes the measurement describing the distribution of intensity within a 3D thresholded region. The region is defined by a threshold on the intensity of the image stack for a given channel.

## Spread
The spread is measured from the [image centered 2nd order moments tensor](https://en.wikipedia.org/wiki/Image_moment) computed in the volume defined by the thresholded area. More specificly, we compute the covariance matrix of the volume. The eigen values of the covariance matrix define the major and minor spread. This is similar to seeing the image as a probablitiy distribution and computing the covariance matrix of the distribution. 

## Entropy
The [entropy](https://en.wikipedia.org/wiki/Entropy_(information_theory)) is computed viewing the image as a probablity density function. The image intensity <img src="https://render.githubusercontent.com/render/math?math=I_p"> is then normalized by the sum over the pixel p in the selection ROI. Finally the entropy is computed as: 
<img src="https://render.githubusercontent.com/render/math?math=H=\sum_{p\in V} \frac{I_p}{\sum_p I_p} \log \frac{I_p}{\sum_p I_p}">
