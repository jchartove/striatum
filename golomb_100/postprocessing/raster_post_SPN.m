function raster_post(directory)
cd(directory);
datadir = [directory, '*data.mat'];
datafiles = dir(datadir);

for file = datafiles'
    filename = strsplit(file.name,'.m');
    filename = filename{1};
    load(file.name);
    
	low_time = 500;
	time_index = time >= low_time;

	FSI_spikes = diff(soma_V(time_index, :) >= 0) == 1;
	D1_spikes = diff(D1_V(time_index, :) >= 0) == 1;
	D2_spikes = diff(D2_V(time_index, :) >= 0) == 1;

	h1 = plot(time(time > low_time)', D2_spikes*diag(1:size(D2_V, 2))', '.', 'Color', [1 .85 0], 'MarkerSize', 10);
	hold on
	h2 = plot(time(time > low_time)', D1_spikes*diag(1:size(D1_V, 2))', '.', 'Color', [.8 .5 .7], 'MarkerSize', 10);
	h3 = plot(time(time > low_time)', FSI_spikes*diag(1:size(soma_V, 2))', '.', 'Color', 'cyan', 'MarkerSize', 10);	%, 'LineWidth', 3)
	xlabel('Time (ms)')
	xlim([500 4000])
	
	mods = simulator_options.modifications;
	val1 = num2str(cell2mat(mods(10,3)')) %9
	val2 = num2str(cell2mat(mods(13,3)')) %11
	val3 = num2str(cell2mat(mods(1,3)'))
	%val3 = num2str(model.fixed_variables.D2_AMPAMSN_AMPA_onset2)
	fulltitle = strcat('strength = ',val1,' freq = ',val2, 'DA lvl = ', val3)
	title(fulltitle)

	saveas(gcf, ['raster_', filename], 'png')
	saveas(gcf, ['raster_', filename], 'fig')
	xlim([1000 1500])
	
	saveas(gcf, ['fig5_zoom_', filename], 'png')
	close all
end
%raster_whoops
