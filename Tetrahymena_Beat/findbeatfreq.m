function [f0, k0] = findbeatfreq(im, s, fps)
n = size(im,1);
m = size(im,2);
%im = im - repmat(mean(im), n, 1);
psd = abs(fft(im));
psd = (psd - repmat(mean(psd), n, 1)) ./ repmat(std(psd), n, 1) ;

g = repmat(1-exp(-([0:n-1]'/s)) , 1, m);
psd = psd .* g;

fmax = max(psd(1:n/2,:),[],2);
fmean = mean(psd(1:n/2,:),2);
f0 = fmax;
k = (0:n/2-1) * fps / n;
k0 = k(find(f0==max(f0),1));
subplot(311);
%imagesc(log(1-min(psd(:))+psd(1:n/2,:))');
imagesc(psd(1:n/2,:)')
subplot(312);
plot(fmax/max(fmax));
hold on
plot(fmean/max(fmean));
hold off
axis tight
title(sprintf('Beat frequency: %.2fHz',k0));
f0 = (f0 - min(f0(:))) ./ (max(f0(:)) - min(f0(:)));
