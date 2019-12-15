%sim_name = 'study_sim1_data';
%DA_level = 'lo';
sim_name = 'study_sim60_data';
DA_level = 'hi';

load(sim_name)
load('FSI_spikes_spectrastats')

figure('Units', 'inches', 'Position', [0 0 6 9.8])

low_time = 500;
time_index = time >= low_time;

%mean_FSI = nanmean(soma_V, 2);
%mean_FSI_detrended = detrend(nanmean(soma_V(time_index, :), 2));
mean_FSI_detrended = detrend(nanmean(soma_soma_somaSomaiSYN_s(time_index, :), 2));

subplot(9, 1, 1)

plot(time(time_index), mean_FSI_detrended, 'LineWidth', 2, 'Color', [0 .8 .8])
axis tight
% set(gca, 'Units', 'inches', 'Position', [0.4, 9.8/7 + .1, 5.6, 9.8/7 - .2])
set(gca, 'Visible', 'off')
xlabel('Simulated LFP')

[s,w,t] = spectrogram(mean_FSI_detrended, 1200, 1100, [0:100], 10000, 'yaxis');

subplot(9, 1, 2)

imagesc(t,w,abs(s))
axis xy
set(gca, 'FontSize', 12, 'XTick', [], 'YTick', (0:20:100))
ylabel('Freq. (Hz)')
pos = get(gca, 'Position');
pos(2) = pos(2) - pos(4);
pos(4) = 2.2*pos(4);
set(gca, 'Position', pos)

subplot(9, 1, 4)
[FSI_hat, F] = pmtm(mean_FSI_detrended,[],[],10000);
% plot(F, FSI_hat, 'LineWidth', 3, 'Color', [0 .8 .8])

if DA_level == 'lo'
datatable = datatable0;
else
datatable = datatable1;
end

plot(mean(datatable), 'LineWidth', 2, 'Color', [0 .8 .8]);
hold on;
errorghost(datatable,1:151,[0 .8 .8]);
    plot(mean(datatable)+std(datatable),'Color','cyan');
    plot(mean(datatable)-std(datatable),'Color','cyan');
set(gca, 'FontSize', 12, 'YTick', [], 'box', 'off')
ylabel('Spectral power')
yyaxis right;
datatable_indiv = indiv(soma_V);
%datatable_indiv = indiv(soma_soma_somaSomaiSYN_s);
plot(mean(datatable_indiv), '--','LineWidth', 2, 'Color', [0 .8 .8]);
yyaxis right;
errorghost(datatable_indiv,1:151,[0 .8 .8]);
    plot(mean(datatable_indiv)+std(datatable_indiv),'Color','cyan');
    plot(mean(datatable_indiv)-std(datatable_indiv),'Color','cyan');
axis('tight');
xlim([0 100])
xlabel('Freq. (Hz)')
set(gca, 'FontSize', 12, 'YTick', [], 'box', 'off')
pos = get(gca, 'Position');
pos(2) = pos(2) - pos(4);
pos(4) = 2.2*pos(4);
set(gca, 'Position', pos)

FSI_spikes = diff(soma_V(time_index, :) >= 0) == 1;

subplot(9, 1, 6)

plot(time(time > low_time)'/1000, FSI_spikes*diag(1:size(soma_V, 2))', '.', 'Color', [0 .8 .8], 'MarkerSize', 10);
ylim([1 size(soma_V, 2)] + .5)
xlim([1 4])
set(gca, 'FontSize', 12, 'YTick', [], 'box', 'off')
pos = get(gca, 'Position');
pos(2) = pos(2) - pos(4);
pos(4) = 2.2*pos(4);
set(gca, 'Position', pos)

subplot(9, 1, 8)

plot(time(time > low_time)'/1000, FSI_spikes*diag(1:size(soma_V, 2))', '.', 'Color', [0 .8 .8], 'MarkerSize', 10);
ylim([1 size(soma_V, 2)] + .5)
xlim([1 1.25])
xlabel('Time (s)')
set(gca, 'FontSize', 12, 'YTick', [], 'box', 'off')
pos = get(gca, 'Position');
pos(2) = pos(2) - pos(4);
pos(4) = 2.2*pos(4);
set(gca, 'Position', pos)

saveas(gcf, ['fig4_', sim_name(1:(end - length('.mat')))])

saveas(gcf, ['fig4_', sim_name(1:(end - length('.mat')))], 'eps')