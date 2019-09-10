figure
subplot1 = subplot(4,2,1);
plot(time,sim9_trace,'k');
ylabel('Voltage (mV)');

subplot2 = subplot(4,2,[2,4]);
p1 = plot(sim9_spectrum, 'LineWidth',2,'DisplayName','I_{app} = 8');
hold on;
p2 = plot(sim21_spectrum, 'LineWidth',2,'DisplayName','I_{app} = 20');
ylabel('Spectral Power (a.u.)');
legend1 = legend([p1 p2]);
set(legend1,'FontSize',16,'AutoUpdate','off','Location','best');
xlim([0 100]);

subplot3 = subplot(4,2,3);
plot(time,sim21_trace,'k');

subplot5 = subplot(4,2,5);
plot(time,sim46_trace,'k');

subplot6 = subplot(4,2,[6,8]);
p3 = plot(mean(datatable_0pt5), 'LineWidth',2,'DisplayName','\lambda = 500');
hold on;
p4 = plot(mean(datatable_7), 'LineWidth',2,'DisplayName','\lambda = 7000');
legend2 = legend([p3 p4]);
set(legend2,'FontSize',16,'AutoUpdate','off','Location','best');
plot(mean(datatable_0pt5)-std(datatable_0pt5),'b');
plot(mean(datatable_0pt5)+std(datatable_0pt5),'b');
plot(mean(datatable_7)-std(datatable_7),'r');
plot(mean(datatable_7)+std(datatable_7),'r');
xlim([0 100]);
xlabel('Frequency (Hz)')

subplot7 = subplot(4,2,7);
plot(time,sim66_trace,'k');
xlabel('Time (ms)')
