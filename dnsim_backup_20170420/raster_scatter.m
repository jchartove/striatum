filenew = strcat(filename, '_FSI')
v_new = soma_v(T_start:T_total+1,:);
[avgfr,spike_pairs, FSI_spike_indicator] = generate_spikes(soma_v, v_new, filenew, time, T_start, dt, numcells);
filenew = strcat(filename, '_D1')
v_new = D1_V(T_start:T_total+1,:);
[avgfr,spike_pairs, D1_spike_indicator] = generate_spikes(D1_V, v_new, filenew, time, T_start, dt, numcells);
%filenew = strcat(filename, '_D2')
%v_new = D2_V(T_start:T_total+1,:);
[avgfr,spike_pairs, D2_spike_indicator] = generate_spikes(D2_V, v_new, filenew, time, T_start, dt, numcells);
[FSI_cell,FSI_time] = find(FSI_spike_indicator)
[D1_cell,D1_time] = find(D1_spike_indicator)
%[D2_cell,D2_time] = find(D2_spike_indicator)
close all

handle5 = figure('units','normalized','outerposition',[0 0 1 1]);
scatter(FSI_time,FSI_cell, 'filled')
hold on
scatter(D1_time,D1_cell,'r', 'filled')
%scatter(D2_time,D2_cell,'g', 'filled')
legend('FSI spikes','D1 SPN spikes')%,'D2 SPN spikes')
xlabel('Time');
ylabel('Cell ID');
imgtitle = strcat(filename,'_scatter.png')
title(imgtitle);
saveas(handle5, imgtitle, 'png');

xlim([(T_start*dt+10000)/10 ((T_start*dt)+30000)/10]);
imgtitle = strcat(filename,'_scatter_zoom.png')
title(imgtitle);
saveas(handle5, imgtitle, 'png');