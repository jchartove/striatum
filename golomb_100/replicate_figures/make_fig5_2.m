%sim_name = 'study_sim1_data';
%DA_level = 'lo';
sim_name = 'study_sim60_data';
DA_level = 'hi';

load(sim_name)

figure('Units', 'inches', 'Position', [0 0 6 9.8])

low_time = 500;
time_index = time >= low_time;

mean_D1 = nanmean(D1_V, 2);
mean_D2 = nanmean(D2_V, 2);

%mean_D1 = nanmean(D1_D1_gabaRecInputMSN_s, 2);
%mean_D2 = nanmean(D2_D2_gabaRecInputMSN_s, 2);

% if strcmp(sim_name, 'study_sim1_data.mat')
%     
%     smoothing = normpdf(-250:.1:250, 0, 100)';
%     smoothing = smoothing/sum(smoothing);
%     
%     mean_D1 = conv(mean_D1, smoothing, 'same');
%     mean_D2 = conv(mean_D2, smoothing, 'same');
%     
% end

D1_spikes = diff(D1_V(time_index, :) >= 0) == 1;
D2_spikes = diff(D2_V(time_index, :) >= 0) == 1;

subplot(9, 1, 8)

h1 = plot(time(time > low_time)', D2_spikes*diag(1:size(D1_V, 2))', '.', 'Color', [1 .85 0], 'MarkerSize', 15);
hold on
h2 = plot(time(time > low_time)', D1_spikes*diag(1:size(D2_V, 2))', '.', 'Color', [.8 .5 .7], 'MarkerSize', 10); %, 'LineWidth', 3)
ylim([1 size(D1_V, 2)] + .5)
set(gca, 'YTick', [], 'FontSize', 12, 'box', 'off') % 'Visible', 'off')
xlabel('Time (ms)')
xlim([1000 4000])
pos = get(gca, 'Position');
pos(2) = pos(2) - pos(4);
pos(4) = 2*pos(4);
set(gca, 'Position', pos)

mean_D2_detrended = detrend(mean_D2(time_index, :));
mean_D1_detrended = detrend(mean_D1(time_index, :));

subplot(9, 1, 1)

[ax, h1, h2] = plotyy(time(time_index), mean_D2_detrended, time(time_index), mean_D1_detrended);
set(h2, 'LineWidth', 2, 'Color', [.8 .5 .7])
set(h1, 'LineWidth', 3, 'Color', [1 .85 0])
legend('D2 SPN', 'D1 SPN')
axis(ax, 'tight')
set(ax, 'Visible', 'off')
% pos = get(gca, 'Position');
% pos(2) = pos(2) + .2*pos(4);
% set(gca, 'Position', pos)

mean_D1_spikes = nanmean(D1_spikes, 2);
[s,w,t] = spectrogram(mean_D1_detrended, 1200, 1100, [0:100], 10000, 'yaxis');
% [s,w,t] = spectrogram(mean_D1_spikes, 1200,1100, [0:100], 10000, 'yaxis');

subplot(9, 1, 2)

imagesc(t,w,abs(s))
axis xy
set(gca, 'FontSize', 12, 'XTick', [], 'YTick', (0:20:100))
ylabel('Freq. (Hz)')
pos = get(gca, 'Position');
pos(2) = pos(2) - pos(4);
pos(4) = 2.2*pos(4);
set(gca, 'Position', pos)
% colorbar

mean_D2_spikes = nanmean(D2_spikes, 2);
[s,w,t] = spectrogram(mean_D2_detrended, 1200, 1100, [0:100], 10000, 'yaxis');
% [s,w,t] = spectrogram(mean_D2_spikes,1200,1100,[0:100],10000,'yaxis');

subplot(9, 1, 4)

imagesc(t,w,abs(s))
axis xy
set(gca, 'FontSize', 12, 'XTick', [], 'YTick', (0:20:100))
ylabel('Freq. (Hz)')
% xlabel('Time (s)')
pos = get(gca, 'Position');
pos(2) = pos(2) - pos(4);
pos(4) = 2.2*pos(4);
set(gca, 'Position', pos)
% colorbar

subplot(9, 1, 6)

%% 

% [D1_hat, F] = pmtm(mean_D1_detrended,[],[],10000);
% plot(F, D1_hat, 'LineWidth', 3, 'Color', [.8 .5 .7]);
% hold on;
% [D2_hat, F] = pmtm(mean_D2_detrended,[],[],10000);
% plot(F, D2_hat, 'LineWidth', 3, 'Color', [1 .85 0]);
%% 

load('D1_spectrastats')
if DA_level == 'lo'
datatable = datatable0;
else
datatable = datatable1;
end
plot(mean(datatable), 'LineWidth', 2, 'Color', [.8 .5 .7]);
hold on;
errorghost(datatable,1:151, [.8 .5 .7]);
    plot(mean(datatable)+std(datatable),'Color',[.8 .5 .7]);
    plot(mean(datatable)-std(datatable),'Color',[.8 .5 .7]);
axis('tight');

load('D2_spectrastats')
if DA_level == 'lo'
datatable = datatable0;
else
datatable = datatable1;
end
plot(mean(datatable), 'LineWidth', 2, 'Color', [1 .85 0]);
    plot(mean(datatable)+std(datatable),'Color',[1 .85 0]);
    plot(mean(datatable)-std(datatable),'Color',[1 .85 0]);
errorghost(datatable,1:151, [1 .85 0]);
%% 

% load('spn_stats.mat');
% 
% DA_level = 'lo';
% 
% if DA_level == 'lo'
% plot(mean(dataD1spikeslow), 'LineWidth', 2, 'Color', [.8 .5 .7]);
% hold on;
% errorghost(dataD1spikeslow,1:151, [.8 .5 .7]);
% else
% plot(mean(dataD1spikeshigh), 'LineWidth', 2, 'Color', [.8 .5 .7]);
% hold on;
% errorghost(dataD1spikeshigh,1:151, [.8 .5 .7]);
% end
% 
 axis('tight');
% 
% if DA_level == 'lo'
% plot(mean(dataD2spikeslow), 'LineWidth', 2, 'Color', [1 .85 0]);
% errorghost(dataD2spikeslow,1:151, [1 .85 0]);
% else
% plot(mean(dataD2spikeshigh), 'LineWidth', 2, 'Color', [1 .85 0]);
% errorghost(dataD2spikeshigh,1:151, [1 .85 0]);
% end
% 
xlim([0 100])
xlabel('Freq. (Hz)')
set(gca, 'FontSize', 12, 'YTick', [], 'box', 'off')
pos = get(gca, 'Position');
pos(2) = pos(2) - pos(4);
pos(4) = 2.2*pos(4);
set(gca, 'Position', pos)
% 
saveas(gcf, ['fig5_', sim_name(1:(end - length('.mat')))]) 
saveas(gcf, ['fig5_', sim_name(1:(end - length('.mat')))], 'eps')
%% 

% figure
% 
% [LFP_hat, F] = pmtm(LFP_trimmed,[],[],10000);
% plot(F, LFP_hat, 'LineWidth', 3, 'Color', 'k')
% xlim([0 100])
% box off
% set(gca, 'FontSize', 16)
% xlabel('Freq. (Hz)')
% 
% saveas(gcf, ['LFP_pmtm_', sim_name(1:(end - length('.mat')))])
% 
% saveas(gcf, ['LFP_pmtm_', sim_name(1:(end - length('.mat')))], 'eps')