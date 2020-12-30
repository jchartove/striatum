function thedecimator(directory)
	cd(directory);
	datadir = [directory, '*.mat'];
	datafiles = dir(datadir);
	
	for file = datafiles'
		filename = strsplit(file.name,'.m');
		filename = filename{1}
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
		time = zeros(1,new_T);
		for j = 1:new_T 
		time(j) = (j-1)*dt;
		end
		clearvars sim_data
		
		save(filename)
		clearvars
		close all
		
	end
end