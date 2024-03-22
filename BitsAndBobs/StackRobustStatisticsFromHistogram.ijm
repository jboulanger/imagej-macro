
function numberOfPixels() {
	Stack.getDimensions(width, height, channels, slices, frames);
	return width * height * channels * slices * frames;
}

function computeHistogramBins(nbins) {
	/*Compute the histogram bins*/
	Stack.getStatistics(count, mean, min, max, stdDev);
	values = newArray(nbins);
	for (i = 0; i < nbins; i++) {
		values[i] = min + (max - min) * i / (nbins);
	}
	return values;
}

function stackHistogram(bins) {
	/* compute the historgram of the all stack */
	Stack.getDimensions(width, height, channels, slices, frames);
	hist = newArray(bins.length);
	for (n = 1; n <= slices; n++) {
		Stack.setSlice(n);
		getHistogram(values, counts, bins.length, bins[0], bins[bins.length-1]);
		for (i = 0; i < nbins; i++) {
			hist[i] = hist[i] + counts[i];
		}
	}	
	return hist;
}

function medianFromHistogram(values, counts) {
	n = 0;
	for (i = 0; i < counts.length; i++) {
		n += counts[i];
	}
	print(n);
	count = 0;
	for (i = 0; i < counts.length; i++) {
		count += counts[i];
		if (count > n/2) {
			med = values[i];
			break;
		}
	}
	return med;
}

function madFromHistogram(values, counts, med) {
	values1 = newArray(values.length);
	mini = 0;
	maxi = maxOf(abs(values[0] - med), abs(values[1] - med));
	n = 0;
	for (i = 0; i < values1.length; i++) {
		values1[i] = mini + (maxi-mini) / nbins * i;
		n += counts[i];
	}
	//  recompute histogram
	counts1 = newArray(nbins);
	for (i = 0; i < values.length; i++) {
		val = abs(values[i] - med);
		idx = (val - mini) / (maxi-mini) * nbins;
		counts1[idx] = counts1[idx] + counts[i];
	}
	//
	count = 0;
	for (i = 0; i < counts1.length; i++) {
		count += counts1[i];
		if (count >= n / 2) {
			mad = values1[i];
			break;
		}
	}
	return 1.4838*mad;
}

function getStackRobustStatictics() {
	/* compute robust statistics median and mad from stack histograms */
	nbins = 10000;	
	values = computeHistogramBins(nbins);
	counts = stackHistogram(values);	
	med = medianFromHistogram(values, counts);
	mad = madFromHistogram(values, counts, med);		
	return newArray(med, mad);
}

print("\\Clear");
run("Close All");
newImage("Untitled", "32-bit black", 500, 500, 10);
run("Add Specified Noise...", "stack standard=10");
stats = getStackRobustStatictics();
Array.print(stats);