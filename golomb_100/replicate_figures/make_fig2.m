%function make_fig2_new()
%CREATEFIGURE(X1, YMatrix1, X2, Y2, X3, Y3)

X1 = iapp;
X2 = taub;
X3 = noise;
YMatrix1(:,1) = hz_no_gd;
YMatrix1(:,2) = hz_gd;
Y2 = burst_taub;
Y3 = burst_noise;

% Create figure
figure
% Create subplot
subplot1 = subplot(1,3,1);
hold(subplot1,'on');

% Create multiple lines using matrix input to plot
plot1 = plot(X1,YMatrix1,'LineWidth',3);
set(plot1(1),'DisplayName','gd = 0');
set(plot1(2),'DisplayName','gd = 6');

% Create ylabel
ylabel('Firing rate within burst (Hz)');

% Create xlabel
xlabel('I_{app} (mA/cm^2)');
%axis(subplot1,'tight');

% Set the remaining axes properties
set(subplot1,'FontSize',16);
legend('Location','southeast');
% Create subplot
subplot2 = subplot(1,3,2);
hold(subplot2,'on');

% Create multiple lines using matrix input to plot
plot(X2,Y2,'LineWidth',3,'DisplayName','\theta freq.');
ylabel('Burst frequency (Hz)');

% Create xlabel
xlabel('tau_{D} (ms)');
axis(subplot2,'tight');
yyaxis right
plot(X2, hz_taub,'LineWidth',3,'DisplayName','\gamma freq.');
ylabel('Firing rate within burst (hz)');
ylim([0 90]);
set(subplot2,'FontSize',16);
legend('Location','southeast');

% Create multiple lines using matrix input to plot
subplot3 = subplot(1,3,3);
hold(subplot3,'on');
plot(X3,mean(Y3),'LineWidth',3,'DisplayName','\theta freq.');
errorghost(Y3,X3,'b');

% Create xlabel
xlabel('Poisson noise rate (MHz)');
ylabel('Burst frequency (Hz)');
yyaxis right
plot(X3,mean(theta_power_noise),'LineWidth',3,'DisplayName','\theta power');
errorghost(theta_power_noise,X3,'r');
legend('Location','southeast');
ylabel('Low frequency power');
axis(subplot3,'tight');
xlim([0 7]);
set(subplot3,'FontSize',16);

% Set the remaining axes properties

%legend1 = legend(subplot1,'show');
%set(legend1,'FontSize',16,'AutoUpdate','off','Location','best');
