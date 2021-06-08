% analyse kymographs
% 
% A csv file list the filenames and the estimated frame rates
% produce by the imageJ macro Beat_Cilia.ijm
%
% read the csv file to get file names and fps
clear all
% load the csv files
% find all csv file
cname = {'WT','MS1','YU3','OAD1C11','Unknown'};% condition names
cpatt = {'WT','M81','YU3','OAD1C11'};% condition pattern 
listcsv = dir('./*/*.csv');
n = 1;
for i = 1:numel(listcsv)
    csvin = readtable(fullfile(listcsv(i).folder, listcsv(i).name));
    for j=1:size(csvin,1)
        % get filename/folder
        str = split(listcsv(i).folder,'/');
        folder{n} = str(end);
        file{n} = csvin.File(j);
        fps(n) = csvin.FrameRate_Hz_(j);
        % get condition
        condidx(n) = numel(cname);
        for k=1:numel(cpatt)
            if ~isempty(strfind(csvin.File{j}, cpatt{k}))
                condidx(n) = k;break;
            end
        end
        % get frequency
        kname = fullfile(listcsv(i).folder, strrep(csvin.File{j},'.nd2','-kymo.tif'));
        im = double(imread(kname));
        [psd(n,:), bps(n)] = findbeatfreq(im, 4, fps(n));
        drawnow;
        n = n + 1;
    end
end
% save results in a csv file
condition = cname(condidx);
result = table(folder',file',condition',fps',bps','VariableName',{'folder','file','condition','fps','bps'});
writetable(result,'result.csv');
%% boxplot
figure(1); clf
plot(condidx+.05*randn(size(condidx)), bps,'.','color',[0.1 0.2 0.8 0.25],'markersize',20);
hold on;
boxplot(bps,condidx);
hold off
xticklabels(cname)
title('Cilia beat frequency')
xlabel('Condition');
ylabel('Beat frequency [Hz]')
axis square
saveas(gcf,'figa-boxplot.pdf')
%% Insight plot
figure(2), clf
for k=1:numel(cname)-1
    subplot(numel(cname)-1,1,k);
    idx = find(condidx == k);
    imagesc(psd(idx,:));
    hold on;
    plot(bps(idx)./fps(idx)*1024, 1:numel(idx),'rs')
    title(cname{k})
end
saveas(gcf,'figb-fft-by-condition.pdf')
%%
figure(3)
hist(fps);
title('Distribution of imaging frame rate');
xlabel('Frame rate [Hz]')
ylabel('Count')
axis square
saveas(gcf,'figc-fps-histogram.pdf')