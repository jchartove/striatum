function generate_img(directory,decim,datatype_range)
	cd(directory);
	datadir = [directory, '*data.mat'];
	datafiles = dir(datadir);
    temp_flag = decim; %kludgey
	
	txtfile = strcat(directory,'.csv');
    txtfile = strrep(txtfile,'/','-')
	formatSpec = '%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s \r\n';
    fileID = fopen(txtfile,'at+');
    tempID5 = fileID; %also kludgey
    fprintf(fileID,formatSpec, ...
    'Filename, Average firing rate, Spike pairs, Total power, Delta, Theta, Alpha, Beta, Low gamma, High gamma, HFO, Low freq peak, Beta peak, Low gamma peak, High gamma peak, HFO peak, Gamma peak, High peak, Checksum, \r\n');
	   
    
	for file = datafiles'
		filename = strsplit(file.name,'.m');
		filename = filename{1};

		%%%%%%%%%decimator code
        if decim
            load(file.name, 'sim_data'); 
            
            [T_total, numcells] = size(sim_data.soma_v);
			numMSNs = size(sim_data.MSN_V,2);
            T_total = T_total-1;
            new_T = T_total/10 + 1;

            soma_v = zeros(new_T,numcells);
            soma_soma_golomb_K_a = zeros(new_T,numcells);
            soma_soma_golomb_K_b = zeros(new_T,numcells);
            soma_soma_golomb_Kdr_n = zeros(new_T,numcells);
            soma_soma_golomb_Na_h = zeros(new_T,numcells);
            soma_soma_soma_soma_iSYN_s = zeros(new_T,numcells);
			
            dend_v = zeros(new_T,numcells);
            dend_dend_golomb_K_a = zeros(new_T,numcells);
            dend_dend_golomb_K_b = zeros(new_T,numcells);
            dend_dend_golomb_Kdr_n = zeros(new_T,numcells);
            dend_dend_golomb_Na_h = zeros(new_T,numcells);
			
 			MSN_V = zeros(new_T,numcells);
            MSN_naCurrentMSN_m = zeros(new_T,numcells);
            MSN_naCurrentMSN_h = zeros(new_T,numcells);
            MSN_kCurrentMSN_m = zeros(new_T,numcells);
 			MSN_mCurrentMSN_m = zeros(new_T,numcells);
 			soma_MSN_soma_MSN_iSYN_s = zeros(new_T,numcells);
			MSN_MSN_gabaRecInputMSN_s = zeros(new_T,numcells);

            for m = 1:numcells
                soma_v(:,m) = decimate(sim_data.soma_v(:,m),10);
                soma_soma_golomb_K_a(:,m) = decimate(sim_data.soma_soma_golomb_K_a(:,m),10); 
                soma_soma_golomb_K_b(:,m) = decimate(sim_data.soma_soma_golomb_K_b(:,m),10);
                soma_soma_golomb_Kdr_n(:,m) = decimate(sim_data.soma_soma_golomb_Kdr_n(:,m),10);
                soma_soma_golomb_Na_h(:,m) = decimate(sim_data.soma_soma_golomb_Na_h(:,m),10);
                soma_soma_soma_soma_iSYN_s(:,m) = decimate(sim_data.soma_soma_soma_soma_iSYN_s(:,m),10);
				
                dend_v(:,m) = decimate(sim_data.dend_v(:,m),10);
                dend_dend_golomb_K_a(:,m) = decimate(sim_data.dend_dend_golomb_K_a(:,m),10);
                dend_dend_golomb_K_b(:,m) = decimate(sim_data.dend_dend_golomb_K_b(:,m),10);
                dend_dend_golomb_Kdr_n(:,m) = decimate(sim_data.dend_dend_golomb_Kdr_n(:,m),10);
                dend_dend_golomb_Na_h(:,m) = decimate(sim_data.dend_dend_golomb_Na_h(:,m),10);
				
			end
			
			for p = 1:numMSNs	
				MSN_V(:,p) = decimate(sim_data.MSN_V(:,p),10);
 				MSN_naCurrentMSN_m(:,p) = decimate(sim_data.MSN_naCurrentMSN_m(:,p),10);
				MSN_naCurrentMSN_h(:,p) = decimate(sim_data.MSN_naCurrentMSN_h(:,p),10);
				MSN_kCurrentMSN_m(:,p) = decimate(sim_data.MSN_kCurrentMSN_m(:,p),10);
				MSN_mCurrentMSN_m(:,p) = decimate(sim_data.MSN_mCurrentMSN_m(:,p),10);
				soma_MSN_soma_MSN_iSYN_s(:,p) = decimate(sim_data.soma_MSN_soma_MSN_iSYN_s(:,p),10);
				MSN_MSN_gabaRecInputMSN_s(:,p) = decimate(sim_data.MSN_MSN_gabaRecInputMSN_s(:,p),10);
            end

            labels = sim_data.labels;
            params = sim_data.params;

            dt = .1;
            clearvars sim_data
			temp_flag = 0;

            save(filename)
        else
            load(file.name);
			temp_flag = 0;
            decim = temp_flag;
        end
		
		%%%%%%%%image generation
		time = zeros(1,new_T);
        for j = 1:new_T 
            time(j) = (j-1)*dt;
        end
		T_total = new_T-1;
        T_start = T_total*0.25; 
		
		for datatype = datatype_range
		
        if datatype == 1
			data = soma_v;
			filenew = strcat(filename, '_FSI')
 		elseif datatype == 2
 			data = MSN_V;
 			filenew = strcat(filename, '_MSN')
		elseif datatype == 3
			data = soma_MSN_soma_MSN_iSYN_s;
			filenew = strcat(filename, '_FSIsyn')
 		elseif datatype == 4
 			data = MSN_mCurrentMSN_m;
 			filenew = strcat(filename, '_mcurr')
 		elseif datatype == 5
 			datat = MSN_MSN_gabaRecInputMSN_s;
 			filenew = strcat(filename, '_MSNsyn')
		end
			
		v_new = data(T_start:T_total+1,:);
		[avgfr,spike_pairs, spike_indicator] = generate_spikes(data, v_new, filenew, time, T_start, dt, numcells);
        fileID = tempID5;
		generate_spec(directory, avgfr, spike_pairs, v_new, filenew, time, dt, numcells, fileID, formatSpec)
		generate_spec(directory, avgfr, spike_pairs, sum(spike_indicator), strcat(filenew, '_spikes'), time, dt, 1, fileID, formatSpec)
			
		%%%%%%%%%%%%%%%%%%%%% gating variables
        handle4 = figure;
		if datatype == 1
			plot(time(T_start+1:end),spike_indicator(1,:),time(T_start+1:end),soma_soma_golomb_Na_h(T_start+1:end,1), ... 
            time(T_start+1:end),soma_soma_golomb_Kdr_n(T_start+1:end,1), time(T_start+1:end),soma_soma_golomb_K_a(T_start+1:end,1), ...
            time(T_start+1:end),soma_soma_golomb_K_b(T_start+1:end,1), time(T_start+1:end), soma_soma_soma_soma_iSYN_s(T_start+1:end,1));
            legend('Spikes','Sodium activation','Potassium activation','Potassium 2 activation', 'Potassium 2 inactivation','IPSC')
 		elseif datatype == 2
             plot(time(T_start+1:end),spike_indicator(1,:),time(T_start+1:end),MSN_naCurrentMSN_h(T_start+1:end,1), ... 
             time(T_start+1:end),MSN_kCurrentMSN_m(T_start+1:end,1), time(T_start+1:end),MSN_naCurrentMSN_m(T_start+1:end,1), ...
             time(T_start+1:end),MSN_mCurrentMSN_m(T_start+1:end,1), time(T_start+1:end), soma_MSN_soma_MSN_iSYN_s(T_start+1:end,1), ...
             time(T_start+1:end),MSN_MSN_gabaRecInputMSN_s(T_start+1:end,1));
             legend('Spikes','Sodium activation','Potassium activation','Sodium inactivation', 'M current activation','FSI to MSN IPSC',...
             'MSN to MSN IPSC')
 		elseif datatype == 4
 			plot(time(T_start+1:end),MSN_mCurrentMSN_m(T_start+1:end,1))
 			legend('M current activation')
		end
        
		xlabel('Time');

        imgtitle = strcat(filenew,'ions.png')
        title(imgtitle);
        saveas(handle4, imgtitle, 'png');

        xlim([T_start*dt+100 (T_start*dt)+200]);
        imgtitle = strcat(filenew,'ions_zoom.png')
        title(imgtitle);
        saveas(handle4, imgtitle, 'png');
		
		close all
		end

        close all
	end
	fclose('all');
end