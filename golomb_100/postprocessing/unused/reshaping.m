fr2 = reshape(fr,[32,8]);

normalized = zeros(size(power));
normalized(:,1) = power(:,1)./(total*3);
normalized(:,2) = power(:,2)./(total*4);
normalized(:,3) = power(:,3)./(total*5);
normalized(:,4) = power(:,4)./(total*23);
normalized(:,5) = power(:,5)./(total*30);
normalized(:,6) = power(:,6)./(total*35);
normalized(:,7) = power(:,7)./(total*50);

normalized2 = zeros(size(power));
normalized2(:,1) = power(:,1)./(total*sqrt(3));
normalized2(:,2) = power(:,2)./(total*sqrt(4));
normalized2(:,3) = power(:,3)./(total*sqrt(5));
normalized2(:,4) = power(:,4)./(total*sqrt(23));
normalized2(:,5) = power(:,5)./(total*sqrt(30));
normalized2(:,6) = power(:,6)./(total*sqrt(35));
normalized2(:,7) = power(:,7)./(total*sqrt(50));

normalized3 = zeros(size(power));
normalized3(:,1) = power(:,1)./(total);
normalized3(:,2) = power(:,2)./(total);
normalized3(:,3) = power(:,3)./(total);
normalized3(:,4) = power(:,4)./(total);
normalized3(:,5) = power(:,5)./(total);
normalized3(:,6) = power(:,6)./(total);
normalized3(:,7) = power(:,7)./(total);

narrow = zeros(size(power));
narrow(:,1) = power(:,1)./3;
narrow(:,2) = power(:,2)./4;
narrow(:,3) = power(:,3)./5;
narrow(:,4) = power(:,4)./23;
narrow(:,5) = power(:,5)./30;
narrow(:,6) = power(:,6)./35;
narrow(:,7) = power(:,7)./50;

% normpeaks = zeros(size(peaks));
% normpeaks(:,1) = peaks(:,1);
% normpeaks(:,2) = peaks(:,2) - 13;
% normpeaks(:,3) = peaks(:,3) - 36;
% normpeaks(:,4) = peaks(:,4) - 66;
% normpeaks(:,5) = peaks(:,5) - 101;
% normpeaks(find(normpeaks < 0)) = 0;

delta = reshape(power(:,1),size(fr2));
theta = reshape(power(:,2),size(fr2));
alpha = reshape(power(:,3),size(fr2));
beta = reshape(power(:,4),size(fr2));
lowgamma = reshape(power(:,5),size(fr2));
highgamma = reshape(power(:,6),size(fr2));
gamma = lowgamma + highgamma;
HFO = reshape(power(:,7),size(fr2));

lowpeak = reshape(peaks(:,1),size(fr2));
bpeak = reshape(peaks(:,2),size(fr2));
lgpeak = reshape(peaks(:,3),size(fr2));
hgpeak = reshape(peaks(:,4),size(fr2));
hfopeak = reshape(peaks(:,5),size(fr2));
gpeak = reshape(peaks(:,6),size(fr2));
hipeak = reshape(peaks(:,7),size(fr2));

total = reshape(total,[32,8]);

delta2 = delta./total;
theta2 = theta./total;
alpha2 = alpha./total;
beta2 = beta./total;
lgamma2 = lowgamma./total;
hgamma2 = highgamma./total;
gamma2 = gamma./total;
HFO2 = HFO./total;

fr3 = divvy(fr2,8);
%fr3 = divvy(fr3',4)';
total3 = divvy(total,8);
%total3 = divvy(total3',4)';

delta3 = divvy(delta2,8);
%delta3 = divvy(delta3',4)';
theta3 = divvy(theta2,8);
%theta3 = divvy(theta3',4)';
alpha3 = divvy(alpha2,8);
%alpha3 = divvy(alpha3',4)';
beta3 = divvy(beta2,8);
%beta3 = divvy(beta3',4)';
lgamma3 = divvy(lgamma2,8);
%lgamma3 = divvy(lgamma3',4)';
hgamma3 = divvy(hgamma2,8);
%hgamma3 = divvy(hgamma3',4)';
gamma3 = divvy(gamma2,8);
%gamma3 = divvy(gamma3',4)';
HFO3 = divvy(HFO2,8);
%HFO3 = divvy(HFO3',4)';

lowpeak2 = divvy(lowpeak,8);
%lowpeak2 = divvy(lowpeak2',4)';
bpeak2 = divvy(bpeak,8);
%bpeak2 = divvy(bpeak2',4)';
lgpeak2 = divvy(lgpeak,8);
%lgpeak2 = divvy(lgpeak2',4)';
hgpeak2 = divvy(hgpeak,8);
%hgpeak2 = divvy(hgpeak2',4)';
hfopeak2 = divvy(hfopeak,8);
%hfopeak2 = divvy(hfopeak2',4)';
gpeak2 = divvy(gpeak,8);
%gpeak2 = divvy(gpeak2',4)';
hipeak2 = divvy(hipeak,8);
%hipeak2 = divvy(hipeak2',4)';

figure
imagesc(fr2)
colorbar
title('Firing rate')
xlabel('taub (D-current)')
ylabel('taumult (spiking potassium)')
figure
imagesc(total)
colorbar
title('Total power')
xlabel('taub (D-current)')
ylabel('taumult (spiking potassium)')

figure
imagesc(delta2)
colorbar
title('Delta power')
xlabel('taub (D-current)')
ylabel('taumult (spiking potassium)')
figure
imagesc(theta2)
colorbar
title('Theta power')
xlabel('taub (D-current)')
ylabel('taumult (spiking potassium)')
figure
imagesc(alpha2)
colorbar
title('Alpha power')
xlabel('taub (D-current)')
ylabel('taumult (spiking potassium)')
figure
imagesc(beta2)
colorbar
title('Beta power')
xlabel('taub (D-current)')
ylabel('taumult (spiking potassium)')
figure
imagesc(lgamma2)
colorbar
title('Low gamma power')
xlabel('taub (D-current)')
ylabel('taumult (spiking potassium)')
figure
imagesc(hgamma2)
colorbar
title('High gamma power')
xlabel('taub (D-current)')
ylabel('taumult (spiking potassium)')
figure
imagesc(gamma2)
colorbar
title('Gamma power')
xlabel('taub (D-current)')
ylabel('taumult (spiking potassium)')
figure
imagesc(HFO2)
colorbar
title('HFO power')
xlabel('taub (D-current)')
ylabel('taumult (spiking potassium)')

figure
imagesc(lowpeak)
colorbar
title('Low frequency peak')
xlabel('taub (D-current)')
ylabel('taumult (spiking potassium)')
figure
imagesc(bpeak)
colorbar
title('Beta frequency peak')
xlabel('taub (D-current)')
ylabel('taumult (spiking potassium)')
figure
imagesc(lgpeak)
colorbar
title('Low gamma frequency peak')
xlabel('taub (D-current)')
ylabel('taumult (spiking potassium)')
figure
imagesc(hgpeak)
colorbar
title('High gamma frequency peak')
xlabel('taub (D-current)')
ylabel('taumult (spiking potassium)')
figure
imagesc(hfopeak)
colorbar
title('HFO frequency peak')
xlabel('taub (D-current)')
ylabel('taumult (spiking potassium)')
figure
imagesc(gpeak)
colorbar
title('Gamma frequency peak')
xlabel('taub (D-current)')
ylabel('taumult (spiking potassium)')
figure
imagesc(hipeak)
colorbar
title('High frequency peak')
xlabel('taub (D-current)')
ylabel('taumult (spiking potassium)')

%%%%%%%%%%%%%
figure
imagesc(fr3)
colorbar
title('Firing rate')
xlabel('taub (D-current)')
ylabel('taumult (spiking potassium)')
figure
imagesc(total3)
colorbar
title('Total power')
xlabel('taub (D-current)')
ylabel('taumult (spiking potassium)')

figure
imagesc(delta3)
colorbar
title('Delta power')
xlabel('taub (D-current)')
ylabel('taumult (spiking potassium)')
figure
imagesc(theta3)
colorbar
title('Theta power')
xlabel('taub (D-current)')
ylabel('taumult (spiking potassium)')
figure
imagesc(alpha3)
colorbar
title('Alpha power')
xlabel('taub (D-current)')
ylabel('taumult (spiking potassium)')
figure
imagesc(beta3)
colorbar
title('Beta power')
xlabel('taub (D-current)')
ylabel('taumult (spiking potassium)')
figure
imagesc(lgamma3)
colorbar
title('Low gamma power')
xlabel('taub (D-current)')
ylabel('taumult (spiking potassium)')
figure
imagesc(hgamma3)
colorbar
title('High gamma power')
xlabel('taub (D-current)')
ylabel('taumult (spiking potassium)')
figure
imagesc(gamma3)
colorbar
title('Gamma power')
xlabel('taub (D-current)')
ylabel('taumult (spiking potassium)')
figure
imagesc(HFO3)
colorbar
title('HFO power')
xlabel('taub (D-current)')
ylabel('taumult (spiking potassium)')

figure
imagesc(lowpeak2)
colorbar
title('Low frequency peak')
xlabel('taub (D-current)')
ylabel('taumult (spiking potassium)')
figure
imagesc(bpeak2)
colorbar
title('Beta frequency peak')
xlabel('taub (D-current)')
ylabel('taumult (spiking potassium)')
figure
imagesc(lgpeak2)
colorbar
title('Low gamma frequency peak')
xlabel('taub (D-current)')
ylabel('taumult (spiking potassium)')
figure
imagesc(hgpeak2)
colorbar
title('High gamma frequency peak')
xlabel('taub (D-current)')
ylabel('taumult (spiking potassium)')
figure
imagesc(hfopeak2)
colorbar
title('HFO frequency peak')
xlabel('taub (D-current)')
ylabel('taumult (spiking potassium)')
figure
imagesc(gpeak2)
colorbar
title('Gamma frequency peak')
xlabel('taub (D-current)')
ylabel('taumult (spiking potassium)')
figure
imagesc(hipeak2)
colorbar
title('High frequency peak')
xlabel('taub (D-current)')
ylabel('taumult (spiking potassium)')