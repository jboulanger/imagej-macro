# Parasite coat recovery

This macro quantifies the time evolution of the surface coating by a fluorescent tagged protein from 3D time lapse image sequences.

The first channel labels the overall structure (inside and outside) of the parasite and the second channel 
on the 3D time lapse experiment labels the signal of interest on the surface.

![image](https://github.com/jboulanger/imagej-macro/assets/3415561/5533b6d2-4fee-4361-83d1-c999485a694d)

After segmentation of the frames of the first channel using an iterative median filter and a threshold, the connected componenet at
the center of the first frame is tracked across the frames of the image sequence whilst others are discarded. The mean
intensity in the mask is then measured over time. Finally, a model taking into the photo-bleaching  is fitted to the measurement
$I(t) = (a + b \ \mathrm{erf} (\frac{t-c}{d}) ) \times e^{-t/e}$ with $c$ the time of the midtime point of the event and $d$ its duration.
The time for the intensity to go from 5% to 95% can be computed as $2.326 \times d$.

## Installation & usage
- Install MorpholibJ by enabling the update site Help>upadte>IJPB-plugins
- Open the image and crop the region you want to analyze with the object of interest centered on the first frame.
- Open the macro in the script editor and press 'run'

## Results

- A label image with the mask of the shell
- A table with the time stamps and mean intensity measurement
- A table with the parameters of the model

|Image|Time mid point [min]|Time constant [min]|5%-95% recovery time [min]|Intensity amplitude [au]|Intensity offset [au]|Photo-bleaching [min]|R squared|
|--|--|--|--|--|--|--|--|
|timelapse.tif|9.307838671|1.596980678|3.714577056|22.514305927|133.144549683|19.980892181|0.998063720|
  
- A figure with the intensity profile and the estimated model.

![image](https://github.com/jboulanger/imagej-macro/assets/3415561/888a3341-c66e-4305-98d0-3121577c7576)


## Visual explanation
The 5% to 95% recovery time can be visualized using the following erf function plot:
```
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.special import erf,erfinv
alpha=0.05
x = np.linspace(-5,5,1000)
y = np.array([0.5+0.5*erf(k) for k in x])
seg = (y > (alpha)) & (y < (1-alpha) )
plt.plot(x,y)
plt.plot(x[seg], 0*np.ones(x[seg].shape),)
plt.plot(x,(alpha)*np.ones(x.shape))
plt.plot(x,1-alpha*np.ones(x.shape))
plt.plot(erfinv(2*alpha-1),alpha,'r.')
plt.plot(erfinv(1-2*alpha),1-alpha,'r.')
plt.xlabel('Time')
plt.ylabel('Error function erf')
plt.title('Recovery')
```

![image](https://github.com/jboulanger/imagej-macro/assets/3415561/69625ac3-7a16-47a4-b4d4-4d56052a1258)
