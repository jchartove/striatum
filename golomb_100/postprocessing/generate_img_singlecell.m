function generate_img_singlecell(directory)
cd(directory);
datadir = [directory, '*data.mat'];
datafiles = dir(datadir);

txtfile = strcat(directory,'.csv');
txtfile = strrep(txtfile,'/','-')
formatSpec = '%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s \r\n';
fileID = fopen(txtfile,'at+');
tempID7 = fileID; %also kludgey
fprintf(fileID,formatSpec, ...
    'Filename, Average firing rate, Spike pairs, Total power, Delta, Theta, Alpha, Beta, Low gamma, High gamma, HFO, Low freq peak, Beta peak, Low gamma peak, High gamma peak, HFO peak, Gamma peak, High peak, Min ISI, Max ISI,  Min ISI 2, Max ISI 2, Checksum, \r\n');


for file = datafiles'
    filename = strsplit(file.name,'.m');
    filename = filename{1};
    load(file.name);
    T_total = size(soma_V,1)-1;
    T_start = T_total*0.25;
    numcells = size(soma_V,2);
    
    %%%%%%%%image generation
    time = zeros(1,size(soma_V,1));
    for j = 1:T_total + 1 
        time(j) = (j-1)*10*simulator_options.dt; %factor of 10 for decimation reasons
    end
    
        data = soma_V;
        filenew = strcat(filename, '_FSI');
        V_short = data(T_start:T_total,:);
        T_in_sec = length(V_short)*simulator_options.dt/100;
        [avgfr,spike_pairs, spike_indicator] = generate_spikes(data, V_short, filenew, time, T_start, simulator_options.dt, numcells);
        fileID = tempID7;
        
        handle5 = figure;
		spike_times = find(spike_indicator == 1);
        ISI = diff(spike_times)*simulator_options.dt/100;
		hist(ISI);
		xlabel('Interspike interval');
        imgtitle = strcat(filenew,'isi.png')
        title(imgtitle);
        saveas(handle5, imgtitle, 'png');
        min_ISI = min(ISI);
        max_ISI = max(ISI);
        
        generate_spec(directory, avgfr, min_ISI, max_ISI, spike_pairs, V_short, filenew, time, simulator_options.dt, numcells, tempID7, formatSpec, simulator_options.modifications)
        generate_spec(directory, avgfr, min_ISI, max_ISI, spike_pairs, spike_indicator, strcat(filenew, '_spikes'), time, simulator_options.dt, 1, tempID7, formatSpec, simulator_options.modifications)
        
        %%%%%%%%%%%%%%%%%%%%% gating Variables
        handle4 = figure;
		plot(time(T_start+1:end),soma_somaGolombNa_h(T_start+1:end,1), ...
        time(T_start+1:end),soma_somaGolombKdr_n(T_start+1:end,1), time(T_start+1:end),soma_somaGolombK_a(T_start+1:end,1), ...
        time(T_start+1:end),soma_somaGolombK_b(T_start+1:end,1));
        legend('Sodium activation','Potassium activation','Potassium 2 activation', 'Potassium 2 inactivation')
               
        xlabel('Time');
        
        imgtitle = strcat(filenew,'ions.png')
        title(imgtitle);
        saveas(handle4, imgtitle, 'png');
        
        xlim([T_start T_start+2000]);
        imgtitle = strcat(filenew,'ions_zoom.png')
        title(imgtitle);
        saveas(handle4, imgtitle, 'png');
        
        %saVe(filename)
        close all
    
    close all
end
fclose('all');
end