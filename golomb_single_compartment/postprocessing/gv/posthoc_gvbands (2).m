function posthoc_gvbands(directory)
	load('studyinfo.mat')
	for n = 1:length(studyinfo.simulations)
		studyinfo.simulations(n).result_functions = {@gvfr};
		studyinfo.simulations(n).result_files = {['study_sim' num2str(n) '_analysis1_gvfr.mat']};
	end
	save('studyinfo.mat','studyinfo')
	data = dsImport(pwd)
	a = dsImportResults('studyinfo.mat',@gvfr)
	for n = 1:length(data)
		dsAnalyze(data(n),{@gvfr},'result_file',['study_sim' num2str(n) '_analysis1_gvfr']);
	end
gv.Run
end

%seperate function for every variable including normalized
%move some preprocessing to main function
%also functions for MSNs

function outdata = gvfr(data)
	fr = dsCalcFR(data,'bin_size',1, 'bin_shift',1);
	fr_soma = {fr.soma_V_FR(1)};
	
	[T_total, numcells] = size(data.soma_V);
	numD1s = size(sim_data.D1_V,2);
	numD2s = size(sim_data.D2_V,2);
    T_total = T_total-1;
    new_T = T_total/10 + 1;
	time = zeros(1,new_T);
    for j = 1:new_T 
        time(j) = (j-1)*dt;
    end
	T_total = size(data.soma_v,1)-1;
    T_start = T_total*0.25; 
	v_new = data.soma_V(T_start:T_total+1,:);

	if numcells > 1
        lfp = mean(v_new');
    else
        lfp = v_new';
    end
    %%%%%%%%%%%%%%%%%%%% spectra
    m = mean(lfp);
    signal = lfp - m; %zero-center
    signal = double(detrend(signal));
	dt = 0.1;
	[y,f] =  pmtm(signal,[],[0:150],1000/dt);
	
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
		
	outdata = {thp};
end