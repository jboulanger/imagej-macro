/*
 * Apply a gamma correction to the current look up table
 */
gamma = getNumber("gamma", 1);

getLut(reds, greens, blues);
for (i = 0; i < reds.length; i++) {
	reds[i] = 255*pow(reds[i]/255,gamma);
	greens[i] = 255*pow(greens[i]/255,gamma);
	blues[i] = 255*pow(blues[i]/255,gamma);
}

setLut(reds, greens, blues);
resetMinAndMax();