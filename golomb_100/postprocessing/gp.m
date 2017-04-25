function fr(directory)
	dt = .01;
	
	cd(directory);
	datadir = [directory, '*.mat'];
	datafiles = dir(datadir);
	txtfile = strcat(directory,'.csv');
	txtfile = strrep(txtfile,'/','-')
	formatSpec = '%s %s %s %s \r\n';
	fileID = fopen(txtfile,'a+');
	fprintf(fileID,formatSpec,'Filename, Gamma power, Gamma peak, raw Gamma power\r\n');
	for file = datafiles'
		filename = strsplit(file.name,'.m');
		filename = filename{1};
		
		load(file.name, 'sim_data'); 
		
		lfp = mean(sim_data.soma_v');
		
		m = mean(lfp);
		signal = lfp - m; %zero-center
		signal = detrend(signal);
		[y] = power_spectrum(signal',sim_data.time,1,1);
		[gp,gi] = max(y(50:100)); %gamma spectral peak
		gi = gi + 50; %giving actual frequency
		gpraw = gp;
		gp = gp/mean(y(10:100)); %normalized. 10-5: took out the lowest frequencies because they're weird
		
		output = {strcat(filename,',') strcat(num2str(gp),',') strcat(num2str(gi),',') strcat(num2str(gpraw),',')}
		fprintf(fileID,formatSpec,output{1,:});
		
        close all
		clearvars sim_data
	end
	fclose('all')
end