function generate_img(directory,decim)
	cd(directory);
	datadir = [directory, '*.mat'];
	datafiles = dir(datadir);
    temp_flag = decim; %kludgey
    
	txtfile = strcat(directory,'.csv');
	txtfile = strrep(txtfile,'/','-')
	formatSpec = '%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s \r\n';
	fileID = fopen(txtfile,'a+');
    tempID = fileID; %also kludgey
	fprintf(fileID,formatSpec, ...
    'Filename, Average firing rate, Spike pairs, Total power, Delta, Theta, Alpha, Beta, Low gamma, High gamma, HFO, Low freq peak, Beta peak, Low gamma peak, High gamma peak, HFO peak, Gamma peak, High peak, \r\n');
    
	for file = datafiles'
		filename = strsplit(file.name,'.m');
		filename = filename{1};

		%%%%%%%%%decimator code
        if decim
            load(file.name, 'sim_data'); 
            
            [T_total, numcells] = size(sim_data.soma_v);
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

            labels = sim_data.labels;
            params = sim_data.params;

            dt = .1;
            clearvars sim_data

            save(filename)
        else
            load(file.name);
            decim = temp_flag;
        end
		
		%%%%%%%%image generation

        time = zeros(1,new_T);
        for j = 1:new_T 
            time(j) = (j-1)*dt;
        end
		T_total = new_T-1;
        T_new = T_total/2; %T_new, not to be confused with new_T. wow gj self
        
        v_new = soma_v(T_new:T_total+1,:); %remove transients
        
        if numcells > 1
            lfpold = mean(soma_v');
            lfp = mean(v_new');
        else
            lfpold = soma_v';
            lfp = v_new';
        end
		
		spike_indicator = zeros(numcells,T_new+1); 
		synch_indicator = zeros(numcells, numcells, T_new);
		
		spikes = zeros(1,numcells);
            
        for t = 1:T_new
            spike_indicator(:,t) = (v_new(t,:)<0) & (v_new(t+1,:) >= 0);
            s = (v_new(t,:)<0) & (v_new(t+1,:) >= 0);
            spikes = spikes + s;
        end
			
        avgfr = mean(spikes)/(T_new/(1000/dt));
		
        spike_pairs = 0;
%         if numcells<=10 % this calculation is so expensive i'm sorry
% 		
%             %rectangle for convolutions
%             rect = 3*ones(1,5/dt);
%             rect = rect(rect > eps);
%             rect = [zeros(1,length(rect)) rect];
% 
%             wide_spikes = zeros(numcells,T_new+1);
%             for c = 1:numcells
%                 wide_spikes(c,:) = conv(spike_indicator(c,:),rect,'same');
%             end
% 
%             for b = 1:numcells
%                 for c = 1:numcells
%                     if b ~= c
%                         wide_sum = wide_spikes(b,:) + wide_spikes(c,:);
%                         foo = wide_sum > 5;
%                         wide_synch = diff(foo, [], 2);
%                         synch_indicator(b,c,:) = wide_synch == 1;
%                         spike_pairs = spike_pairs + sum(synch_indicator(b,c,:));
%                         c
%                     end
%                 end
%             end
%         end

		%%%%%%%%%%%%%%%%%%%% spike plots
		handle0 = figure;
		hist(spikes/(T_new/(1000/dt)))
		xlabel('Firing rate');
		imgtitle = strcat(filename,'hist.png')
        title(imgtitle);
		saveas(handle0, imgtitle, 'png');
		
		handle1 = figure;
		plot(time,soma_v);
		hold on;
		plot(time, lfpold,'k','LineWidth',2);
		hold off;
		xlabel('Time');
		ylabel('Voltage');
		imgtitle = strcat(filename,'spikes.png')
        title(imgtitle);
		saveas(handle1, imgtitle, 'png');
		
		xlim([T_new*dt (T_new*dt)+100]);
		imgtitle = strcat(filename,'spikes_zoom.png')
        title(imgtitle);
		saveas(handle1, imgtitle, 'png');
		
		%%%%%%%%%%%%%%%%%%%% rasters
		handle2 = figure;
		imagesc(soma_v');
		colorbar;
		xlabel('Time');
		ylabel('Cell number');
		imgtitle = strcat(filename,'raster.png')
        title(imgtitle);
		saveas(handle2, imgtitle, 'png');
		
		xlim([T_new T_new+(100/dt)]);
		imgtitle = strcat(filename,'raster_zoom.png')
        title(imgtitle);
		saveas(handle2, imgtitle, 'png');
		
		%%%%%%%%%%%%%%%%%%%% spectrograms
		handle3 = figure;
		
		m = mean(lfp);
		signal = lfp - m; %zero-center
		signal = detrend(signal);
		[y] = power_spectrum(signal',time,0,0);
		totalp = sum(y(1:150)); %total power. below: eeg bands
        dp = sum(y(1:3));
        thp = sum(y(4:7));
        ap = sum(y(8:12));
        
        %for the broader peaks, also find peak location
		[~,lowpeak] = max(y(1:12));
        bp = sum(y(13:35));
        [~,bpeak] = max(y(13:35));
        bpeak = bpeak + 12;
		gplow = sum(y(36:65));
        [~,glopeak] = max(y(36:65));
        glopeak = glopeak + 35;
		gphigh = sum(y(66:100));
        [~,ghipeak] = max(y(66:100));
        ghipeak = ghipeak + 65;
        hfop = sum(y(101:150));
        [~,hfopeak] = max(y(101:150));
        hfopeak = hfopeak + 100;
		[~,gpeak] = max(y(36:100));
        gpeak = gpeak + 35;
		[~,hipeak] = max(y(66:150));
        hipeak = hipeak + 65;
		
		spectitle = strcat(filename,'spectrum.mat')
		save(spectitle,'y');
		imgtitle = strcat(filename,'spectrum.png')
        title(imgtitle);
		saveas(handle3, imgtitle, 'png');
		
		%%%%%%%%%%%%%%%%%%%%% gating variables
		handle4 = figure;
		plot(time(T_new+1:end),spike_indicator(1,:),time(T_new+1:end),soma_soma_golomb_Na_h(T_new+1:end,1), ... 
			time(T_new+1:end),soma_soma_golomb_Kdr_n(T_new+1:end,1), time(T_new+1:end),soma_soma_golomb_K_a(T_new+1:end,1), ...
			time(T_new+1:end),soma_soma_golomb_K_b(T_new+1:end,1), time(T_new+1:end), soma_soma_soma_soma_iSYN_s(T_new+1:end,1));
		legend('Spikes','Sodium activation','Potassium activation','Potassium 2 activation', 'Potassium 2 inactivation','IPSC')
		xlabel('Time');
		
		imgtitle = strcat(filename,'ions.png')
        title(imgtitle);
		saveas(handle4, imgtitle, 'png');
		
		xlim([T_new*dt (T_new*dt)+100]);
		imgtitle = strcat(filename,'ions_zoom.png')
        title(imgtitle);
		saveas(handle4, imgtitle, 'png');
 		
        fileID = tempID;
		output = {strcat(filename,',') strcat(num2str(avgfr),',') strcat(num2str(spike_pairs),',') ... 
            strcat(num2str(totalp),',') strcat(num2str(dp),',') strcat(num2str(thp),',') strcat(num2str(ap),',') ...
            strcat(num2str(bp),',') strcat(num2str(gplow),',') strcat(num2str(gphigh),',') strcat(num2str(hfop),',')  ...
            strcat(num2str(lowpeak),',') strcat(num2str(bpeak),',') strcat(num2str(glopeak),',') strcat(num2str(ghipeak),',') ...
			strcat(num2str(hfopeak),',') strcat(num2str(gpeak),',') strcat(num2str(hipeak),',')};
		fprintf(fileID,formatSpec,output{1,:});
		
        close all
	end
	fclose('all');
end