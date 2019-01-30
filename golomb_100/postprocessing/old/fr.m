function fr(directory)
	dt = .01;
	
	cd(directory);
	datadir = [directory, '*.mat'];
	datafiles = dir(datadir);
	txtfile = strcat(directory,'.csv');
	txtfile = strrep(txtfile,'/','-')
	formatSpec = '%s %s\r\n';
	fileID = fopen(txtfile,'a+');
	fprintf(fileID,formatSpec,'Filename, Average firing rate\r\n');
	for file = datafiles'
		filename = strsplit(file.name,'.m');
		filename = filename{1};
		
		load(file.name, 'sim_data'); 
		[T, numcells] = size(sim_data.soma_v);
		T = T-1;
		
		spikes = zeros(1,numcells);

		for t = 1:T
			s = (sim_data.soma_v(t,:)<0) & (sim_data.soma_v(t+1,:) >= 0);
			spikes = spikes + s;
		end
		avgfr = mean(spikes)/(T/(1000/dt))

		output = {strcat(filename,',') strcat(num2str(avgfr),',')};
		fprintf(fileID,formatSpec,output{1,:});
		
        close all
		clearvars sim_data
	end
	fclose('all')
end