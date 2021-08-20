% Fit a Gaussian function to an airy pattern
%
% 
clear all

lambda = 600;
NA = 1.45;


%% What is the relationshop between the std of a gaussian and the rayleight criterion
%
% Fit a Gaussian to an Airy function
disp('*** Rayleigh and Gaussian STD ***');
r = linspace(-5,5,10000);
% define functions
airy = @(x) (2*real(besselj(1,x))./x).^2;
gauss = @(p,x) (exp(-.5*(x/p(1)).^2));
cosexp = @(p,x) 0.5*(1+cos(2*pi*x/p(1))).*exp(-.5*(x/(p(2))).^2);
if ~exist('OCTAVE_VERSION','var')
     nonlin_curvefit = @(a,b,c,d) lsqcurvefit(a,b,c,d,[],[],optimset('display','off'));     
   
end
localmin = @(x,y) x( y<circshift(y,-1) & y<circshift(y,1) );

% Estimate the standard devition of the matching Gaussian by fitting a Gaussian 
% to an Airy pattern
sigma = nonlin_curvefit(gauss, 1, r, airy(r));

% Estimate the standard devition of the matching Sinc exp by fitting it 
% to an Airy pattern
a = nonlin_curvefit(cosexp, [8;pi/sqrt(2)], r, airy(r));

% Find the value of the 1st zeros of the Airy function
z = localmin(r,airy(r));
z = min(z(z>0));
fwhm = (max(r(airy(r)>=.5))-min(r(airy(r)>.5)));
fprintf('Airy FWHM %.2fnm (%.4f lambda/NA)\n', fwhm * lambda / (2*pi*NA), round(fwhm / (2*pi)*1000)/1000);
fprintf('Rayleigh %.2fnm (%.4f lambda/NA)\n',z * lambda / (2*pi*NA), round(z / (2*pi)*1000)/1000);
fprintf('Gaussian std %.2fnm (%f lambda/NA)\n', sigma * lambda / (2*pi*NA), sigma / (2*pi));
fprintf('cos exp [%.2fnm %.2fnm]\n', a * lambda / (2*pi*NA));

fprintf('ratio FWHM/std %f\n', fwhm / sigma);
fprintf('ratio rayleigh/std %f\n', z / sigma);
% Plot the 3 functions
figure(1);
plot(r,airy(r),'r');
hold on;
plot(r,gauss(sigma,r),'b');
hold on
plot(r,cosexp(a,r),'g');
plot([-fwhm/2 fwhm/2],[0.5,0.5],'k');
plot([0 z],[0.3,0.3],'k');
hold off
legend('airy','gaussian','sinc*gaussian')
title('1d simulation');

%% How does a Gaussian fit is influence by the size of the bead
%
figure(2);
disp('*** Bead size and Gaussian PSF fit ***');

% Convolve a 2D bead with a 2D airy pattern
% fit a gaussian function and compare to a prediction
px = 4;
N = 512;
d = 300;
vals = ((-N/2-1):N/2)*px+.5;
[x,y] = meshgrid(vals,vals);
D = double(sqrt(x.^2+y.^2) < d/2);
H = airy(sqrt(x.^2+y.^2)/lambda * (pi*NA));
%I = fftshift(real(ifft2(fft2(D).*fft2(H))));
I = conv2(H,D,'same');
I = I ./ max(I(:));
sigmab = nonlin_curvefit(gauss, 500, vals, I(256,:));
sigmap = sqrt((sigma * lambda / (pi*NA))^2 + (d/4).^2);
%a = sigma * lambda / (pi*NA);
%b = d/4;
%sigmap = sqrt(a^2*b^2/(a^2-b^2));
fprintf('Gaussian fit std %f\n', sigmab);
fprintf('Gaussian est std %f\n', sigmap);
fprintf('Relative error %.2f%%\n', (sigmap-sigmab)/(sigmab) * 100);
plot(x(256,:),sqrt(I(256,:)),'r');
hold on;
plot(x(256,:),sqrt(H(256,:)),'g');
plot(x(256,:),sqrt(D(256,:)),'m');
plot(vals,sqrt(gauss(sigmab,vals)),'b');
hold off
legend('intensity','psf','bead','gaussian')
title('2d simulation with bead')



%% What percentage of modulation corresponds to the Rayleigh criterion
%
% use the MTF to compute the modulation amplitude corresponding to the Rayleigh
% and FWHM criterion
% 
% Use coltman formula for the CTF
disp('*** MTF, Rayleight and FWHM ***');
figure(3);
otf = @(v) 2 / pi * real(acos(abs(v)) - abs(v) .* sqrt(max(0,1-v.^2)));
ctf = @(v) 4 / pi * (otf(v)-otf(3*v)/3+otf(5*v)/5-otf(7*v)/7+otf(9*v)/9-otf(11*v)/11+otf(13*v)/13-otf(15*v)/15+otf(17*v)/17);
k = linspace(0,1);
plot(k,otf(k),'linewidth',2);
hold on
plot(k,ctf(k),'m','linewidth',2);
plot(k,1-k/.8,'g','linewidth',2);
plot(k,1-k,'g','linewidth',2);
k1 = 0.5/0.52;
plot([0 k1 k1],[otf(k1), otf(k1) 0],'r','linewidth',2);
k2 = 0.5/0.61;
plot([0 k2 k2],[otf(k2), otf(k2) 0],'r','linewidth',2);
hold off;
grid on
axis([0 1 0 1])
title('MTF, CTF, Rayleigh and FWHM')
legend('MTF','CTF','1-k')
fprintf('modulation amplitude at the fwhm crit.     : %f%%\n',otf(k1)*100)
fprintf('modulation amplitude at the rayleigh crit. : %f%%\n',otf(k2)*100)


% 
px = 100;
N = 512;
d = 100;
vals = ((-N/2-1):N/2)*px+.5;
[x,y] = meshgrid(vals,vals);
H = airy(sqrt(x.^2+y.^2)/lambda * (pi*NA));
H = H/sum(H(:));

clear c
for f=1:256
%D = double(mod(x/p,1)<0.5); 
D = 0.5*sin(2*pi*f*x/px);
I = fftshift(real(ifft2(fft2(D).*fft2(H))));
%I = conv2(D,H,'same');
%I = I(N/4:3*N/4,N/4:3*N/4);
%plot(D(N/2,:)); hold on; plot(I(N/2,:),'r');hold off; drawnow
c(f) = max(I(:))-min(I(:)) / (max(I(:))+min(I(:)));
end
plot(c)

