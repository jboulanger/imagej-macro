/*
 * T test 
 * 
 * Compare with matlab implementation
 * dat = importdata('data.csv'); [h,p,ci,stats]=ttest2(dat.data(:,1),dat.data(:,2),'left')
 * Compare with R
 * df<-read.csv('data.csv')
 * t.test(df$x1,df$x2,"less")
 */
 
 
// create 2 random vectors
x1 = generateSamples(10,0,1);
x2 = generateSamples(10,0,1);

// display the data in a table
Array.show("data",x1,x2);
saveAs("Results", "/home/jeromeb/data.csv");
ttest(x1,x2,"less");


function generateSamples(n,m,s) {
	/* Generate a Normal vector 
	 * input
	 * n: length
	 * m: mean
	 * s: standard deviation
	 */
	x = newArray(n);
	for (i = 0; i < n; i++) {
		x[i] = m+s*random("gaussian");
	}
	return x;
}

///////////////// T Test /////////////////////

function ttest(x1,x2,tail) {
	/* independent two-sample t-test
	 *  input:
	 *  x1, x2 : arrays
	 *  tail: tail: two sided, less, greater
	 *  
	 *  TODO implement tail (for now it is only less)
	 */
	Array.getStatistics(x1, min, max, m1, s1);
	n1 = x1.length;
	Array.getStatistics(x2, min, max, m2, s2);
	n2 = x2.length;
	if (n1==n2 && (s1<2*s2||s2<2*s1)) { // equal sample size, equal variance
		s = sqrt((s1*s1+s2*s2)/2);
		t = (m1 - m2) / (s * sqrt(2 / n1));
		df = 2 * n1 - 1;
		case = "equal size, equal variance";
	} else { 
		if (s1<2*s2||s2<2*s1) {
			s = sqrt(((n1-1)*s1*s1+(n2-1)*s2*s2) / (n1+n2-2));
			t = (m1 - m2) / (s*sqrt(1/n1+1/n2));
			df = n1+n2-2;
			case = "non-equal size, equal variance";
		} else {
			s = sqrt(s1*s1/n1+s2*s2/n2);
			t = (m1 - m2) / s;
			df = pow(s1*s1/n1+s2*s2/n2,2) / ( pow(s1*s1/n1,2)/(n1-1) +  pow(s2*s2/n2,2)/(n2-1));
			case = "non-equal size, non-equal variance";
		}
	}
	p = student_t_cdf(t, df);
	
	print("hypt:"+tail+ " case:" + case + " t:"+t + ", df:"+df + " p:"+p);
	return p;
}

function lgamma(z) {
	/* Log Gamma function
	 * Gergő Nemes approximation of the Log Gamma function 
	 */
	return 0.5 * ( log(2*PI) - log(z)) + z * (log(z + 1 / (12 * z - 1/(10*z))) - 1);
}

function incbeta(a,b,x) {
	/* Incomplete Beta function
	 * https://en.wikipedia.org/wiki/Beta_function#Incomplete_beta_function
	 * adatpted from https://github.com/codeplea/incbeta/blob/master/incbeta.c
	 */
	 
	if (x < 0.0 || x > 1.0) return 1.0/0.0;
	if (x > (a+1.0)/(a+b+2.0)) {
        return (1.0-incbeta(b,a,1.0-x)); /*Use the fact that beta is symmetrical.*/
    }
    lbeta_ab = lgamma(a)+lgamma(b)-lgamma(a+b);
    front = exp(log(x)*a+log(1.0-x)*b-lbeta_ab) / a;
    f = 1.0;
    c = 1.0; 
    d = 0.0;
    
    for (i = 0; i <= 200; ++i) {
    	 m = i/2;
    	 if (i == 0) {
            numerator = 1.0; /*First numerator is 1.0.*/
        } else if (i % 2 == 0) {
            numerator = (m*(b-m)*x)/((a+2.0*m-1.0)*(a+2.0*m)); /*Even term.*/
        } else {
            numerator = -((a+m)*(a+b+m)*x)/((a+2.0*m)*(a+2.0*m+1)); /*Odd term.*/
        }
        d = 1.0 + numerator * d;
        if (abs(d) < 1e-8) {d = 1e-8;}
        d = 1.0 / d;

        c = 1.0 + numerator / c;
        if (abs(c) < 1e-8) {c = 1e-8};

        cd = c*d;
        f *= cd;

        /*Check for stop.*/
        if (abs(1.0-cd) < 1e-30) {
            return front * (f-1.0);
        }    
    }
    return 1.0/0.0;
}

function student_t_cdf(t, v) {
    /* The cumulative distribution function (CDF) for Student's t distribution
     *  
    */
    x = (t + sqrt(t * t + v)) / (2.0 * sqrt(t * t + v));
    prob = incbeta(v/2.0, v/2.0, x);
    return prob;
}