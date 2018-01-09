function sta(directory)
cd(directory);
datadir = [directory, '*data.mat'];
datafiles = dir(datadir);
windowsize = 2000;

for file = datafiles'
    filename = strsplit(file.name,'.m');
    filename = filename{1};
    load(file.name, 'simulator_options','D1_V','D1_D1_gabaRecInputMSN_s','D1_soma_somaMSNiSYN_s','D1_mCurrentMSN_m','soma_V');
    T_total = size(D1_V,1)-1;
    T_start = T_total*0.25;
    new_T = T_total + 1; %removed a factor of 10. this is all very silly at this point
    numcells = size(D1_V,2);
    
    time = zeros(1,size(D1_V));
    for j = 1:new_T
        time(j) = (j-1)*simulator_options.dt;
    end
    filenew = strcat(filename, '_D1STA')
    V_new = D1_V(T_start:T_total+1,:);
    T_new = length(time)-T_start-1; %this is dumb but saves me time rewriting
    spike_indicator = zeros(numcells,T_new+1);
    
    spikes = zeros(1,numcells);
    
    for t = 1:T_new
        spike_indicator(:,t) = (V_new(t,:)<0) & (V_new(t+1,:) >= 0);
        s = (V_new(t,:)<0) & (V_new(t+1,:) >= 0);
        spikes = spikes + s;
    end
    
    SPN_in = D1_D1_gabaRecInputMSN_s(5001:end,:)';
    FSI_in = D1_soma_somaMSNiSYN_s(5001:end,:)';
    m_in = D1_mCurrentMSN_m(5001:end,:)';
    V = D1_V(5001:end,:)';
    SPN_pre = zeros(numcells,windowsize);
    SPN_ntwk = zeros(numcells,windowsize);
    FSI_pre = zeros(numcells,windowsize);
    FSI_ntwk = zeros(numcells,windowsize);
    m_pre = zeros(numcells,windowsize);
    m_ntwk = zeros(numcells,windowsize);
    V_pre = zeros(numcells,windowsize);
    V_ntwk = zeros(numcells,windowsize);
    for i = 1:numcells
        [SPN_pre(i,:), ~] = evTrigAvg(spike_indicator(i,:),SPN_in(i,:), windowsize);
        [FSI_pre(i,:), ~] = evTrigAvg(spike_indicator(i,:),FSI_in(i,:), windowsize);
        [m_pre(i,:), ~] = evTrigAvg(spike_indicator(i,:),m_in(i,:), windowsize);
        [V_pre(i,:), ~] = evTrigAvg(spike_indicator(i,:),V(i,:), windowsize);
        [SPN_ntwk(i,:), ~] = evTrigAvg(spike_indicator(i,:),mean(SPN_in), windowsize);
        [FSI_ntwk(i,:), ~] = evTrigAvg(spike_indicator(i,:),mean(FSI_in), windowsize);
        [m_ntwk(i,:), ~] = evTrigAvg(spike_indicator(i,:),mean(m_in), windowsize);
        [V_ntwk(i,:), ~] = evTrigAvg(spike_indicator(i,:),mean(V), windowsize);
    end
    close all
    save(filenew);
    %%%%%%%%%%%%%%%%%%%% plots
%     handle1 = figure;
%     plot((-200:0.1:-0.1),nanmean(m_pre));
%     hold on;
%     plot((-200:0.1:-0.1),nanmean(m_ntwk),'k','LineWidth',2);
%     hold off;
%     xlabel('Time before D1 spike (ms)');
%     ylabel('M current');
%     imgtitle = strcat(filenew,'m.png')
%     title(imgtitle);
%     saveas(handle1, imgtitle, 'png');
%     
%     handle2 = figure;
%     plot((-200:0.1:-0.1),nanmean(SPN_pre));
%     hold on;
%     plot((-200:0.1:-0.1),nanmean(SPN_ntwk),'k','LineWidth',2);
%     hold off;
%     xlabel('Time before D1 spike (ms)');
%     ylabel('SPN input');
%     imgtitle = strcat(filenew,'spn.png')
%     title(imgtitle);
%     saveas(handle2, imgtitle, 'png');
%     
%     handle3 = figure;
%     plot((-200:0.1:-0.1),nanmean(FSI_pre));
%     hold on;
%     plot((-200:0.1:-0.1),nanmean(FSI_ntwk),'k','LineWidth',2);
%     hold off;
%     xlabel('Time before D1 spike (ms)');
%     ylabel('FSI input');
%     imgtitle = strcat(filenew,'fsi.png')
%     title(imgtitle);
%     saveas(handle3, imgtitle, 'png');
    
    handle4 = figure;
    plot((-200:0.1:-0.1),nanmean(V_pre));
    hold on;
    plot((-200:0.1:-0.1),nanmean(V_ntwk),'k','LineWidth',2);
    hold off;
    xlabel('Time before D1 spike (ms)');
    ylabel('D1 voltage');
    imgtitle = strcat(filenew,'v.png')
    title(imgtitle);
    saveas(handle4, imgtitle, 'png');
end

fclose('all');
end