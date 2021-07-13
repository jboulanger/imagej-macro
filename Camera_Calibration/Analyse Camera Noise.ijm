//@String(Label="Name",value="camera") name
//@Integer(Label="Radius",value=4) radius
//@Boolean(Label="Robust",value=true) robust
//@Float(Label="DC offset (<0:estimate)",value=-1) dcmean
//@Integer(Label="Samples",value=10000) n
//@Boolean(Label="Curve in photon",value=true) inphoton
//@Float(Label="Qe",value=0.95) Qe
//@Float(Label="Pixel Size [nm]",value=65) pixelsize
/*
 * Analyse Camera Noise
 * 
 * Acquire a stack of images with different intensity levels 
 * 
 * Ideal sample is a Chroma slide
 * 
 * Retreive the camera gain, dark level offset and standard deviation
 * 
 * The number of photo-electron is (I-offset)/gain
 * The noise variance is gain * I + noise offset
 * The readout noise in electron is (noise std) / gain
 * 
 * Jerome Boulanger 2018
 */
 
macro "Analyse Camera Noise" {
	id=getImageID;			
	E = newArray(n*nSlices);
    V = newArray(n*nSlices);  
    setBatchMode(true);  
    
    computeLocalStatistics(E,V,radius,robust, n);	
			
	// Adjust a 1st order polymomial (a+bx)
    Fit.doFit("Straight Line", E, V);
    gain = Fit.p(1);
    edc = Fit.p(0);    
	selectImage(id);
	
    setSlice(nSlices);	
	run("Select All");

	// if the dark current offset is not provided, estimate it
	if (dcmean < 1) {
		Array.getStatistics(E, dcmean, max, mean, stdDev);
		Array.getStatistics(V, dcstd, max, mean, stdDev);
		dcstd = sqrt(dcstd);
	} else {
		dcstd = gain * dcmean + edc;
		if (dcstd < 0) {
			print("Warning: the value of dcmean (" + dcmean + ") is probably too low");
			print("         It should be at least "+(-edc/gain));
			dcstd = 0;
		} else {
			dcstd = sqrt(gain * dcmean + edc);	
		}
	}
	row = nResults;
	setResult("System",row, name);
	setResult("A/D Gain",row, gain);	
	setResult("Readout [e-]",row, dcstd/gain);
	setResult("DC offset",row, dcmean);
	setResult("DC std",row, dcstd);		
	setResult("Noise offset (calc.)",row, dcstd * dcstd-gain * dcmean);	
	setResult("Noise offset (fit.)",row, edc);	
	setResult("R2",row, Fit.rSquared);
	setResult("Qe",row, Qe);
	setResult("Pixel Size",row, pixelsize);	
	run("Select None");
	if (inphoton) {
		plotNoiseSignatureInPhoton(E,V,gain,dcmean,dcstd,edc,Qe);
	} else {
		plotNoiseSignature(E,V,gain,edc);
	}
	//Plot.create("Curve ", "Intensity", "Noise Variance",E,V);
	setBatchMode(false);
	//Fit.doFit("2nd Degree Polynomial", E, V);
	//Fit.plot();
}

function computeLocalStatistics(E,V,radius,robust,n) {
	id = getImageID();
 	name = getTitle; 	 	
 	// compute local statistics (mean and noise variance) using filters 	
 	if (robust==false) { 		
 		run("Duplicate...", "duplicate title=Mean");
 		run("32-bit");
		run("Mean...", "radius="+radius+" stack");
		selectImage(id);
		run("Duplicate...", "duplicate title=Variance");
		run("32-bit");
		//run("Convolve...", "text1=[0.0 0.22361 0.0\n0.22361 -0.89443 0.22361\n0.0 0.22361 0.0\n] stack");		
		//run("Convolve...", "text1=[0.22361 0 0.22361\n0 -0.89443 0\n0.22361 0 0.22361 \n] stack");
		run("Convolve...", "text1=[0.5 -0.5 0\n-0.5 0.5 0\n0 0 0\n] stack");
		run("Variance...", "radius="+radius+ " stack");
 	} else { 		
 		run("Duplicate...", "duplicate title=Mean");
 		run("32-bit");
 		run("Median...", "radius="+radius+" stack");
		selectImage(id);
		// Median of Absolute Deviation
		run("Duplicate...", "duplicate title=Variance");
		run("32-bit");
		//run("Convolve...", "text1=[0.0 0.22361 0.0\n0.22361 -0.89443 0.22361\n0.0 0.22361 0.0\n] stack");
		//run("Convolve...", "text1=[0.22361 0 0.22361\n0 -0.89443 0\n0.22361 0.0 0.22361\n] stack");
		run("Convolve...", "text1=[0.5 -0.5 0\n-0.5 0.5 0\n0 0 0\n] stack");
		run("Duplicate...", "duplicate title=tmp");
		run("Median...", "radius="+radius+" stack");
		imageCalculator("Subtract", "Variance", "tmp");
		run("Abs","stack");		
		run("Median...", "radius="+radius+" stack");
		run("Multiply...", "value=1.4826 stack");
		run("Square","stack");
		selectWindow("tmp"); close();
 	}	
    
	// sample n measure points in the image (too long otherwise)	  
	x =  newArray(n);    
	y =  newArray(n);  
	W = getWidth();
	H = getHeight();
	for (k = 0; k < nSlices; k++) {		
		//get a set of n point location
		for (i = 0; i < n; i++) {	        
       		x[i] = round(radius + random()*(W-2*radius));
           	y[i] = round(radius + random()*(H-2*radius));            
        }
        
        // measure the intensity in the mean image for these n location
      	selectWindow("Mean");      	
      	//setSlice(k+1);      	
      	setZCoordinate(k);      	
      	for (i = 0; i < n; i++) {  	       
           E[i+k*n] = getPixel(x[i],y[i]);           
      	}
      	
        // measure the intensity in the variance image for these n location   
        selectWindow("Variance");       
        //setSlice(k+1);
        setZCoordinate(k);
        for (i = 0; i < n; i++) {       	                   
            V[i+k*n] = getPixel(x[i],y[i]);     
        }
	}	
}

function plotNoiseSignature(E,V,gain,offset) {
        Array.getStatistics(E,Emin,Emax,a,b);
        Array.getStatistics(V,Vmin,Vmax,a,b);
        //print("Emin="+Emin+", Emax="+Emax);
        // prepare data to draw the model
        samples = 10;
        Ef = newArray(samples);
        Vf = newArray(samples);
        for (i = 0; i < samples; i++){
        	Ef[i] = Emin + i * (Emax-Emin) / (samples-1);
        	Vf[i] = gain * Ef[i] + offset;
        }
        // Create a graph
        Plot.create("Noise Analysis ["+getTitle+"]", "Intensity", "Noise Variance");
        Plot.setLimits(Emin, Emax, Vmin, Vmax);
        Plot.setColor("#5dade2");
        Plot.add("dots", E, V);
        Plot.setColor("#ffaa00");
        Plot.add("line", Ef, Vf);
        Plot.setLegend("Data\nFit");               
}

/*
 * Plot the noise signature expressed in photons
*/
function plotNoiseSignatureInPhoton(E,V,gain,dcmean,dcstd,offset,Qe) {	
	   // convert intensities to photons	   
	   En = newArray(E.length);
       Vn = newArray(V.length);        
	   for (i = 0; i < E.length; i++) {
        	En[i] = (E[i] - dcmean) / gain / Qe;
        	Vn[i] = V[i] / (gain*gain);
        }                  
        Array.getStatistics(En,Emin,Emax,a,b);
        Array.getStatistics(Vn,Vmin,Vmax,a,b);
        // prepare data to draw the model
        samples = 10;
        Ef = newArray(samples);
        Vf = newArray(samples);        
        for (i = 0; i < samples; i++){
        	Ef[i] = Emin + i * (Emax-Emin) / (samples-1);        	
        	Vf[i] = Ef[i] + (dcstd)/(gain*gain);
        }
        // Create a graph       
        Plot.create("Noise Analysis ["+getTitle+"]", "Photon count [e-]", "Noise Variance [e-]");
        Plot.setLimits(Emin, Emax, Vmin, Vmax);
        Plot.setColor("#5dade2");
        Plot.add("dots", En, Vn);
        Plot.setColor("#ffaa00");
        Plot.add("line", Ef, Vf);
        Plot.setLegend("Data\nFit");          
}