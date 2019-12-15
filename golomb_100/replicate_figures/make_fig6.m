load(sim_name)

figure('Units', 'inches', 'Position', [0 0 6 9.8*(5/7)])

low_time = 500;
time_index = time >= low_time;

mean_FSI = nanmean(soma_soma_somaSomaiSYN_s, 2);
mean_D1 = 5*nanmean(D1_D1_gabaRecInputMSN_s, 2)+nanmean(D1_soma_somaMSNiSYN_s, 2);
mean_D2 = 5*nanmean(D2_D2_gabaRecInputMSN_s, 2)+nanmean(D2_soma_somaMSNiSYN_s, 2);

% if strcmp(sim_name, 'study_sim1_data.mat')
%     
%     smoothing = normpdf(-250:.1:250, 0, 100)';
%     smoothing = smoothing/sum(smoothing);
%     
%     mean_D1 = conv(mean_D1, smoothing, 'same');
%     mean_D2 = conv(mean_D2, smoothing, 'same');
%     
% end

mean_D2_detrended = detrend(mean_D2(time_index, :));
mean_D1_detrended = detrend(mean_D1(time_index, :));
mean_FSI_detrended = detrend(mean_FSI(time_index, :));

% subplot(15, 1, 4)
% 
% plot(time(time_index), mean_FSI_detrended, 'LineWidth', 2, 'Color', 'c')
% axis tight
% set(gca, 'Visible', 'off')
% pos = get(gca, 'Position');
% pos(2) = pos(2) - .1*pos(4);
% pos(4) = 1.2*pos(4);
% set(gca, 'Position', pos)
% 
% subplot(15, 1, 5)
% 
% plot(time(time_index), mean_D1_detrended, 'LineWidth', 2, 'Color', [.8 .5 .7])
% axis tight
% set(gca, 'Visible', 'off')
% pos = get(gca, 'Position');
% pos(2) = pos(2) - .1*pos(4);
% pos(4) = 1.2*pos(4);
% set(gca, 'Position', pos)
% 
% subplot(15, 1, 6)
% 
% plot(time(time_index), mean_D2_detrended, 'LineWidth', 2, 'Color', [1 .85 0])
% axis tight
% set(gca, 'Visible', 'off')
% pos = get(gca, 'Position');
% pos(2) = pos(2) - .1*pos(4);
% pos(4) = 1.2*pos(4);
% set(gca, 'Position', pos)

LFP = mean_FSI_detrended + mean_D1_detrended + mean_D2_detrended;
LFP_trimmed = LFP; %(time_index, :);

subplot(4, 5, [1,2,3,4])

plot(time(time_index), LFP_trimmed, 'LineWidth', 2, 'Color', 'k')
axis tight
%set(gca, 'Visible', 'off')
pos = get(gca, 'Position');
pos(2) = pos(2) - .1*pos(4);
pos(4) = 1.2*pos(4);
set(gca, 'Position', pos)

[s,w,t] = spectrogram(LFP_trimmed,1200,1100,[0:100],10000,'yaxis');

subplot(4, 5, [6,7,8,9,11,12,13,14])

imagesc(t,w,abs(s))
axis xy
set(gca, 'FontSize', 12, 'XTick', [], 'YTick', (0:20:100))
ylabel('Freq. (Hz)')
%pos = get(gca, 'Position');
%pos(2) = pos(2) - 1.2*pos(4);
%pos(4) = 2.2*pos(4);
%set(gca, 'Position', pos)

beta = s(w == 18, :);
gamma = s(w == 80, :);

subplot(4, 5, [16,17,18,19])

[ax, h1, h2] = plotyy(t, abs(beta), t, abs(gamma));
axis(ax, 'tight')
set(h1, 'LineWidth', 3)
set(h2, 'LineWidth', 3)
legend({'\beta Power', '\gamma Power'})
%xlim([0 1.5])
%set(gca, 'FontSize', 12, 'XTick', [0 .5 1 1.5])
xlabel('Time (s)')
pos = get(gca, 'Position');
pos(4) = 1.2*pos(4);
set(gca, 'Position', pos)
set(ax, 'box', 'off')
set(ax, 'YTickLabel', [])

%figure('Units', 'inches', 'Position', [0 0 9.8*(3/7) 2])

subplot(4, 5, [10,15])

[LFP_hat, F] = pmtm(LFP_trimmed,[],[],10000);
plot(LFP_hat, F, 'LineWidth', 3, 'Color', 'k')
axis tight
ylim([0 100])
box off
set(gca, 'FontSize', 12)
set(gca, 'visible', 'off')
ylabel('Freq. (Hz)')


saveas(gcf, ['fig6_', sim_name(1:(end - length('.mat')))])

saveas(gcf, ['fig6_', sim_name(1:(end - length('.mat')))], 'eps')

%saveas(gcf, ['LFP_pmtm_', sim_name(1:(end - length('.mat')))])

%saveas(gcf, ['LFP_pmtm_', sim_name(1:(end - length('.mat')))], 'eps')