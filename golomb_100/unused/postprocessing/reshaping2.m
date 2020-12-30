num_rows = 21;
num_cols = 11;
xlabel_str = 'Applied current'; %columns variable
xlims = [0 10];
ylabel_str = 'Noisy input rate'; %rows variable
ylims = [0 4];

fr2 = reshape(fr,[num_rows,num_cols]);
total2 = reshape(total,[num_rows,num_cols]);

normalized = zeros(size(power));
normalized2 = zeros(size(power));
normalized3 = zeros(size(power));
narrow = zeros(size(power));

numbands = 8;
starts = [3,4,5,23,30,35,50];
names = {'Delta', 'Theta', 'Alpha', 'Beta', 'Low_gamma', 'High_gamma', 'HFO', 'Gamma'};
peaknames = {'Low_frequency_peak','Beta_peak','Low_gamma_peak','High_gamma_peak','HFO_peak','Gamma_peak','High_peak'};

figure
titlestr = 'Firing rate.png';
handle1 = imagesc(xlims, ylims, fr2);
colormap('jet');
axis xy
colorbar
title(titlestr)
xlabel(xlabel_str)
ylabel(ylabel_str)
saveas(handle1, titlestr)

figure
titlestr = 'Total power.png';
handle2 = imagesc(xlims, ylims, total2);
colormap('jet');
axis xy
colorbar
title(titlestr)
xlabel(xlabel_str)
ylabel(ylabel_str)
saveas(handle2, titlestr)

close all

for i = 1:numbands
    if i < 7
        bands.(char(names(i))).power = reshape(power(:,i),size(fr2));
        bands.(char(peaknames(i))) = reshape(peaks(:,i),size(fr2));
        
        figure
        titlestr = [char(peaknames(i)) '.png'];
        handle3 = imagesc(xlims, ylims, bands.(char(peaknames(i))));
        colormap('jet');
        axis xy
        colorbar
        title(titlestr)
        xlabel(xlabel_str)
        ylabel(ylabel_str)
        saveas(handle3, titlestr)
        
        normalized(:,i) = power(:,i)./(total*starts(i));
        normalized2(:,i) = power(:,i)./(total*sqrt(starts(i)));
        normalized3(:,i) = power(:,i)./(total);
        narrow(:,i) = power(:,i)./starts(i);
    else 
        bands.(char(names(i))).power = bands.(char(names(5))).power + bands.(char(names(6))).power; %total gamma
    end
    
    bands.(char(names(i))).normalized = bands.(char(names(i))).power./total2; 
    
    figure
	titlestr = ['Raw ' char(names(i)) ' power.png'];
	handle2 = imagesc(xlims, ylims, bands.(char(names(i))).power);
    colormap('jet');
	axis xy
	colorbar
	title(titlestr)
	xlabel(xlabel_str)
	ylabel(ylabel_str)
	saveas(handle2, titlestr)
    
    figure
	titlestr = [char(names(i)) ' power.png'];
	handle1 = imagesc(xlims, ylims, bands.(char(names(i))).normalized);
    colormap('jet');
	axis xy
	colorbar
	title(titlestr)
	xlabel(xlabel_str)
	ylabel(ylabel_str)
	saveas(handle1, titlestr)

	close all
end
	
