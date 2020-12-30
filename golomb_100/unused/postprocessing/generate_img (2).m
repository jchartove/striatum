function generate_img(directory,decim,datatype_range)
	cd(directory);
	datadir = [directory, '*data.mat'];
	datafiles = dir(datadir);
    temp_flag = decim; %kludgey
	
	txtfile = strcat(directory,'.csv');
    txtfile = strrep(txtfile,'/','-')
	formatSpec = '%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s \r\n';
    fileID = fopen(txtfile,'at+');
    tempID7 = fileID; %also kludgey
    fprintf(fileID,formatSpec, ...
    'Filename, Average firing rate, Spike pairs, Total power, Delta, Theta, Alpha, Beta, Low gamma, High gamma, HFO, Low freq peak, Beta peak, Low gamma peak, High gamma peak, HFO peak, Gamma peak, High peak, Checksum, \r\n');
	   
    
	for file = datafiles'
		filename = strsplit(file.name,'.m');
		filename = filename{1};

		%%%%%%%%%decimator code
        if decim
            load(file.name, 'sim_data'); 
            
            [T_total, numcells] = size(sim_data.soma_v);
			numD1s = size(sim_data.D1_V,2);
			numD2s = size(sim_data.D2_V,2);
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
			
 			D1_V = zeros(new_T,numD1s);
            D1_naCurrentMSN_m = zeros(new_T,numD1s);
            D1_naCurrentMSN_h = zeros(new_T,numD1s);
            D1_kCurrentMSN_m = zeros(new_T,numD1s);
 			D1_mCurrentMSN_m = zeros(new_T,numD1s);
 			soma_D1_soma_MSN_iSYN_s = zeros(new_T,numD1s);
			D1_D1_gabaRecInputMSN_s = zeros(new_T,numD1s);
			
 			D2_V = zeros(new_T,numD2s);
            D2_naCurrentMSN_m = zeros(new_T,numD2s);
            D2_naCurrentMSN_h = zeros(new_T,numD2s);
            D2_kCurrentMSN_m = zeros(new_T,numD2s);
 			D2_mCurrentMSN_m = zeros(new_T,numD2s);
 			soma_D2_soma_MSN_iSYN_s = zeros(new_T,numD2s);
			D2_D2_gabaRecInputMSN_s = zeros(new_T,numD2s);

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
			
			for p = 1:numD1s	
				D1_V(:,p) = decimate(sim_data.D1_V(:,p),10);
 				D1_naCurrentMSN_m(:,p) = decimate(sim_data.D1_naCurrentMSN_m(:,p),10);
				D1_naCurrentMSN_h(:,p) = decimate(sim_data.D1_naCurrentMSN_h(:,p),10);
				D1_kCurrentMSN_m(:,p) = decimate(sim_data.D1_kCurrentMSN_m(:,p),10);
				D1_mCurrentMSN_m(:,p) = decimate(sim_data.D1_mCurrentMSN_m(:,p),10);
				soma_D1_soma_MSN_iSYN_s(:,p) = decimate(sim_data.soma_D1_soma_MSN_iSYN_s(:,p),10);
				D1_D1_gabaRecInputMSN_s(:,p) = decimate(sim_data.D1_D1_gabaRecInputMSN_s(:,p),10);
            end
			
			for p = 1:numD2s	
				D2_V(:,p) = decimate(sim_data.D2_V(:,p),10);
 				D2_naCurrentMSN_m(:,p) = decimate(sim_data.D2_naCurrentMSN_m(:,p),10);
				D2_naCurrentMSN_h(:,p) = decimate(sim_data.D2_naCurrentMSN_h(:,p),10);
				D2_kCurrentMSN_m(:,p) = decimate(sim_data.D2_kCurrentMSN_m(:,p),10);
				D2_mCurrentMSN_m(:,p) = decimate(sim_data.D2_mCurrentMSN_m(:,p),10);
				soma_D2_soma_MSN_iSYN_s(:,p) = decimate(sim_data.soma_D2_soma_MSN_iSYN_s(:,p),10);
				D2_D2_gabaRecInputMSN_s(:,p) = decimate(sim_data.D2_D2_gabaRecInputMSN_s(:,p),10);
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
		T_total = size(soma_v,1)-1;
        T_start = T_total*0.25; 
		
		for datatype = datatype_range
		
        if datatype == 1
			data = soma_v;
			filenew = strcat(filename, '_FSI')
 		elseif datatype == 2
 			data = D1_V;
 			filenew = strcat(filename, '_D1')
		elseif datatype == 3
			data = soma_D1_soma_MSN_iSYN_s;
			filenew = strcat(filename, '_FSID1syn')
 		elseif datatype == 4
 			data = D1_mCurrentMSN_m;
 			filenew = strcat(filename, '_D1mcurr')
 		elseif datatype == 5
 			data = D1_D1_gabaRecInputMSN_s;
 			filenew = strcat(filename, '_D1syn')
		elseif datatype == 6
 			data = D2_V;
 			filenew = strcat(filename, '_D2')
		elseif datatype == 7
			data = soma_D2_soma_MSN_iSYN_s;
			filenew = strcat(filename, '_FSID2syn')
 		elseif datatype == 8
 			data = D2_mCurrentMSN_m;
 			filenew = strcat(filename, '_D2mcurr')
 		elseif datatype == 9
 			data = D2_D2_gabaRecInputMSN_s;
 			filenew = strcat(filename, '_D2syn')
		end
			
		v_new = data(T_start:T_total+1,:);
		[avgfr,spike_pairs, spike_indicator] = generate_spikes(data, v_new, filenew, time, T_start, dt, numcells);
        %theres a bug here i haven't fixed due to laziness
        %it should only be numcells in the line above if data is soma_v.
        %should change for d1s and d2s. whatever
        fileID = tempID7;
		generate_spec(directory, avgfr, spike_pairs, v_new, filenew, time, dt, numcells, tempID7, formatSpec)
		generate_spec(directory, avgfr, spike_pairs, sum(spike_indicator), strcat(filenew, '_spikes'), time, dt, 1, tempID7, formatSpec)
			
		%%%%%%%%%%%%%%%%%%%%% gating variables
        handle4 = figure;
		if datatype == 1
			plot(time(T_start+1:end),spike_indicator(1,:),time(T_start+1:end),soma_soma_golomb_Na_h(T_start+1:end,1), ... 
            time(T_start+1:end),soma_soma_golomb_Kdr_n(T_start+1:end,1), time(T_start+1:end),soma_soma_golomb_K_a(T_start+1:end,1), ...
            time(T_start+1:end),soma_soma_golomb_K_b(T_start+1:end,1), time(T_start+1:end), soma_soma_soma_soma_iSYN_s(T_start+1:end,1));
            legend('Spikes','Sodium activation','Potassium activation','Potassium 2 activation', 'Potassium 2 inactivation','IPSC')
 		elseif datatype == 2
             plot(time(T_start+1:end),spike_indicator(1,:),time(T_start+1:end),D1_naCurrentMSN_h(T_start+1:end,1), ... 
             time(T_start+1:end),D1_kCurrentMSN_m(T_start+1:end,1), time(T_start+1:end),D1_naCurrentMSN_m(T_start+1:end,1), ...
             time(T_start+1:end),D1_mCurrentMSN_m(T_start+1:end,1), time(T_start+1:end), soma_D1_soma_MSN_iSYN_s(T_start+1:end,1), ...
             time(T_start+1:end),D1_D1_gabaRecInputMSN_s(T_start+1:end,1));
             legend('Spikes','Sodium activation','Potassium activation','Sodium inactivation', 'M current activation','FSI to D1 IPSC',...
             'D1 to D1 IPSC')
 		elseif datatype == 4
 			plot(time(T_start+1:end),D1_mCurrentMSN_m(T_start+1:end,1))
 			legend('M current activation')
 		elseif datatype == 6
             plot(time(T_start+1:end),spike_indicator(1,:),time(T_start+1:end),D2_naCurrentMSN_h(T_start+1:end,1), ... 
             time(T_start+1:end),D2_kCurrentMSN_m(T_start+1:end,1), time(T_start+1:end),D2_naCurrentMSN_m(T_start+1:end,1), ...
             time(T_start+1:end),D2_mCurrentMSN_m(T_start+1:end,1), time(T_start+1:end), soma_D2_soma_MSN_iSYN_s(T_start+1:end,1), ...
             time(T_start+1:end),D2_D2_gabaRecInputMSN_s(T_start+1:end,1));
             legend('Spikes','Sodium activation','Potassium activation','Sodium inactivation', 'M current activation','FSI to D2 IPSC',...
             'D2 to D2 IPSC')
 		elseif datatype == 8
 			plot(time(T_start+1:end),D2_mCurrentMSN_m(T_start+1:end,1))
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
		
		%save(filename)
		close all
		end

        close all
	end
	fclose('all');
end