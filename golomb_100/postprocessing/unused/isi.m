function isi(directory)
cd(directory);
datadir = [directory, '*data.mat'];
datafiles = dir(datadir);

txtfile = strcat(directory,'.csv');
txtfile = strrep(txtfile,'/','-')
formatSpec = '%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s \r\n';
fileID = fopen(txtfile,'at+');
tempID7 = fileID; %also kludgey
fprintf(fileID,formatSpec, ...
    'Filename, Average firing rate, Spike pairs, Total power, Delta, Theta, Alpha, Beta, Low gamma, High gamma, HFO, Low freq peak, Beta peak, Low gamma peak, High gamma peak, HFO peak, Gamma peak, High peak, Min ISI, Max ISI,  Min ISI 2, Max ISI 2, Checksum, \r\n');


for file = datafiles'
    filename = strsplit(file.name,'.m');
    filename = filename{1};
    load(file.name);
    if exist('FSI_V','var')
        soma_V = FSI_V;
    end
    T_total = size(soma_V,1)-1;
    T_start = T_total*0.25;
    numcells = size(soma_V,2);
    dt = simulator_options.dt;
    mods = simulator_options.modifications;
    
    %%%%%%%%image generation
    time = zeros(1,size(soma_V,1));
    for j = 1:T_total + 1
        time(j) = (j-1)*10*simulator_options.dt; %factor of 10 for decimation reasons
    end
    
    data = soma_V;
    filenew = strcat(filename, '_FSI');
    V_short = data(T_start:T_total,:);
    T_in_sec = length(V_short)*simulator_options.dt/100;
    
    
    T_short = length(V_short)-1;
    if numcells > 1
        lfp = mean(data');
    else
        lfp = data';
    end
    
    spike_indicator = zeros(numcells,T_short);
    spikes = zeros(1,numcells);
    
    for t = 1:T_short
        spike_indicator(:,t) = (V_short(t,:)<0) & (V_short(t+1,:) >= 0);
        s = (V_short(t,:)<0) & (V_short(t+1,:) >= 0);
        spikes = spikes + s;
    end
    
    T_in_sec = (T_short)*dt/100; %this is 100 for decimation reasons
    
    avgfr = mean(spikes)/T_in_sec;
    
    spike_times = find(spike_indicator == 1);
    ISI = diff(spike_times)*dt/100;
    
    min_ISI = min(ISI);
    max_ISI = max(ISI);
    
    m = mean(lfp);
    signal = lfp - m; %zero-center
    signal = double(detrend(signal));
    [y,f] =  pmtm(signal',[],[0:150],10000);
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
    
    checksum = 0;
    
    output = {strcat(filenew,',') strcat(num2str(avgfr),',')...
        strcat(num2str(totalp),',') strcat(num2str(dp),',') strcat(num2str(thp),',') strcat(num2str(ap),',') ...
        strcat(num2str(bp),',') strcat(num2str(gplow),',') strcat(num2str(gphigh),',') strcat(num2str(hfop),',')  ...
        strcat(num2str(lowpeak),',') strcat(num2str(bpeak),',') strcat(num2str(glopeak),',') strcat(num2str(ghipeak),',') ...
        strcat(num2str(hfopeak),',') strcat(num2str(gpeak),',') strcat(num2str(hipeak),',')  ...
        strcat(num2str(min_ISI),',') strcat(num2str(max_ISI),',')  ...
        strcat(num2str(checksum),',') strcat(strjoin(mods(:,1)'),',') strcat(strjoin(mods(:,2)'),',') num2str(cell2mat(mods(:,3)')) }
    fprintf(fileID,formatSpec,output{1,:});
    close all
end
fclose('all');
end