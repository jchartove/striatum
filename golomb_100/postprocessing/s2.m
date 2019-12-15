datadir = pwd;
datafiles = dir(strcat(datadir,'/*data.mat'));
indices = struct;
indices.files = zeros(length(datafiles),1);
indices.onset1 = zeros(length(datafiles),1);
indices.onset2 = zeros(length(datafiles),1);
indices.spikes_in_A = zeros(length(datafiles),35000);
indices.spikes_in_B = zeros(length(datafiles),35000);
indices.spikes_not_in_A = zeros(length(datafiles),35000);
indices.spikes_not_in_B = zeros(length(datafiles),35000);
indices.index_A = zeros(length(datafiles),35000);
indices.index_B = zeros(length(datafiles),35000);
indices.theta_angle = zeros(length(datafiles),35000);
indices.beta_angle = zeros(length(datafiles),35000);
iterator = 0;
%time_index = 500:0.1:4000;

for file = datafiles'
    iterator = iterator + 1

    filenum = strsplit(file.name,{'study_sim','_data'},'CollapseDelimiters',true);
    indices.files(iterator) = str2num(filenum{2});
    load(file.name);
    
    mods = simulator_options.modifications;
	indices.onset1(iterator) = cell2mat(mods(11,3)');
	indices.onset2(iterator) = model.fixed_variables.D1_AMPAMSN_AMPA_onset2;
    
    low_time = 500;
time_index = time >= low_time;
    D1_spikes = diff(D1_V(time_index, :) >= 0) == 1;
    
    %all_spikes = movsum(sum(D1_spikes'),200);
    indices.spikes_in_A(iterator,:) = movsum(sum(D1_spikes(:,51:100)'),200);
    indices.spikes_in_B(iterator,:) = movsum(sum(D1_spikes(:,26:75)'),200);
    indices.spikes_not_in_A(iterator,:) = movsum(sum(D1_spikes(:,1:50)'),200);
    indices.spikes_not_in_B(iterator,:) = movsum(sum(D1_spikes(:,1:25)')+sum(D1_spikes(:,76:100)'),200);
    indices.index_A(iterator,:) = (indices.spikes_in_A(iterator,:)-indices.spikes_not_in_A(iterator,:));%./all_spikes;
    indices.index_B(iterator,:) = (indices.spikes_in_B(iterator,:)-indices.spikes_not_in_B(iterator,:));%./all_spikes;
    
    lfp = mean(D1_V(time_index, :)');
    m = mean(lfp);
    signal = lfp - m; %zero-center
    signal = double(detrend(signal));
%     [LFP_hat, F] = pmtm(signal,[],[],10000);
    
    theta = wavelet_spectrogram(signal, 10000, 3, 2, 1);
    beta = wavelet_spectrogram(signal, 10000, 18, 6, 1);
    
    theta_angle1 = angle(theta);
    beta_angle1 = angle(beta);
    indices.theta_angle(iterator,:) = theta_angle1(1:length(indices.index_A(iterator,:)));
    indices.beta_angle(iterator,:) = beta_angle1(1:length(indices.index_A(iterator,:)));
    %scatter(theta_angle(1:length(index_A)),index_A)
    close all;
end
fclose('all');
T = struct2table(indices); % convert the struct array to a table
sortedT = sortrows(T, 'files'); % sort the table
sortedS = table2struct(sortedT); % change it back to struct array

%% 

theta_angle = {sortedS.theta_angle}.';
index_A = {sortedS.index_A}.';
beta_angle = {sortedS.beta_angle}.';
index_B = {sortedS.index_B}.';
theta_angle = cell2mat(theta_angle);
index_A = cell2mat(index_A);
beta_angle = cell2mat(beta_angle);
index_B = cell2mat(index_B);

spikes_in_A = {sortedS.spikes_in_A}.';
spikes_in_B = {sortedS.spikes_in_B}.';
spikes_in_A = cell2mat(spikes_in_A);
spikes_in_B = cell2mat(spikes_in_B);

spikes_not_in_A = {sortedS.spikes_not_in_A}.';
spikes_not_in_B = {sortedS.spikes_not_in_B}.';
spikes_not_in_A = cell2mat(spikes_not_in_A);
spikes_not_in_B = cell2mat(spikes_not_in_B);

new_index_A = spikes_in_A.*(spikes_in_A > spikes_not_in_A);
new_index_B = spikes_in_B.*(spikes_in_B > spikes_not_in_B);
%% 

%range_no_FS = 1:2:271;
%range_FS = 273:2:543;
%range_no_FS = 1:272;
%range_FS = 273:544;
range_no_FS = 1:130;
range_FS = 131:260;

theta_angle_no_FS = theta_angle(range_no_FS,:);
theta_angle_FS = theta_angle(range_FS,:);
index_A_no_FS = index_A(range_no_FS,:);
index_A_FS = index_A(range_FS,:);

beta_angle_no_FS = beta_angle(range_no_FS,:);
beta_angle_FS = beta_angle(range_FS,:);
index_B_FS = index_B(range_FS,:);
index_B_no_FS = index_B(range_no_FS,:);

new_index_A_no_FS = new_index_A(range_no_FS,:);
new_index_A_FS = new_index_A(range_FS,:);
new_index_B_FS = new_index_B(range_FS,:);
new_index_B_no_FS = new_index_B(range_no_FS,:);

spikes_in_A_no_FS = spikes_in_A(range_no_FS,:);
spikes_in_A_FS = spikes_in_A(range_FS,:);
spikes_in_B_no_FS = spikes_in_B(range_no_FS,:);
spikes_in_B_FS = spikes_in_B(range_FS,:);

theta_angle_no_FS_2 = reshape(theta_angle_no_FS,1,length(range_FS)*35000);
theta_angle_FS_2 = reshape(theta_angle_FS,1,length(range_FS)*35000);
index_A_no_FS_2 = reshape(index_A_no_FS,1,length(range_FS)*35000);
index_A_FS_2 = reshape(index_A_FS,1,length(range_FS)*35000);

beta_angle_no_FS_2 = reshape(beta_angle_no_FS,1,length(range_FS)*35000);
beta_angle_FS_2 = reshape(beta_angle_FS,1,length(range_FS)*35000);
index_B_no_FS_2 = reshape(index_B_no_FS,1,length(range_FS)*35000);
index_B_FS_2 = reshape(index_B_FS,1,length(range_FS)*35000);

new_index_A_no_FS_2 = reshape(new_index_A_no_FS,1,length(range_FS)*35000);
new_index_A_FS_2 = reshape(new_index_A_FS,1,length(range_FS)*35000);
new_index_B_no_FS_2 = reshape(new_index_B_no_FS,1,length(range_FS)*35000);
new_index_B_FS_2 = reshape(new_index_B_FS,1,length(range_FS)*35000);

spikes_in_A_no_FS_2 = reshape(spikes_in_A_no_FS,1,length(range_FS)*35000);
spikes_in_A_FS_2 = reshape(spikes_in_A_FS,1,length(range_FS)*35000);
spikes_in_B_no_FS_2 = reshape(spikes_in_B_no_FS,1,length(range_FS)*35000);
spikes_in_B_FS_2 = reshape(spikes_in_B_FS,1,length(range_FS)*35000);

theta_A_FS = table(theta_angle_FS_2',index_A_FS_2',spikes_in_A_FS_2',new_index_A_FS_2','VariableNames',{'Phase','Index','Spikes','Newindex'});
theta_B_FS = table(theta_angle_FS_2',index_B_FS_2',spikes_in_B_FS_2',new_index_B_FS_2','VariableNames',{'Phase','Index','Spikes','Newindex'});
theta_A_no_FS = table(theta_angle_no_FS_2',index_A_no_FS_2',spikes_in_A_no_FS_2',new_index_A_no_FS_2','VariableNames',{'Phase','Index','Spikes','Newindex'});
theta_B_no_FS = table(theta_angle_no_FS_2',index_B_no_FS_2',spikes_in_B_no_FS_2',new_index_B_no_FS_2','VariableNames',{'Phase','Index','Spikes','Newindex'});
beta_A_FS = table(beta_angle_FS_2',index_A_FS_2',spikes_in_A_FS_2',new_index_A_FS_2','VariableNames',{'Phase','Index','Spikes','Newindex'});
beta_B_FS = table(beta_angle_FS_2',index_B_FS_2',spikes_in_B_FS_2',new_index_B_FS_2','VariableNames',{'Phase','Index','Spikes','Newindex'});
beta_A_no_FS = table(beta_angle_no_FS_2',index_A_no_FS_2',spikes_in_A_no_FS_2',new_index_A_no_FS_2','VariableNames',{'Phase','Index','Spikes','Newindex'});
beta_B_no_FS = table(beta_angle_no_FS_2',index_B_no_FS_2',spikes_in_B_no_FS_2',new_index_B_no_FS_2','VariableNames',{'Phase','Index','Spikes','Newindex'});

s_theta_A_FS = sortrows(theta_A_FS, 'Phase');
s_theta_B_FS = sortrows(theta_B_FS, 'Phase');
s_theta_A_no_FS = sortrows(theta_A_no_FS, 'Phase');
s_theta_B_no_FS = sortrows(theta_B_no_FS, 'Phase');
s_beta_A_FS = sortrows(beta_A_FS, 'Phase');
s_beta_B_FS = sortrows(beta_B_FS, 'Phase');
s_beta_A_no_FS = sortrows(beta_A_no_FS, 'Phase');
s_beta_B_no_FS = sortrows(beta_B_no_FS, 'Phase');

n_theta_A_FS = s_theta_A_FS;
n_theta_B_FS = s_theta_B_FS; 
n_theta_A_no_FS = s_theta_A_no_FS; 
n_theta_B_no_FS = s_theta_B_no_FS; 
n_beta_A_FS = s_beta_A_FS; 
n_beta_B_FS = s_beta_B_FS; 
n_beta_A_no_FS = s_beta_A_no_FS; 
n_beta_B_no_FS = s_beta_B_no_FS; 

toDelete = n_theta_A_FS.Newindex ==0;
n_theta_A_FS(toDelete,:) = [];
toDelete = n_theta_B_FS.Newindex ==0;
n_theta_B_FS(toDelete,:) = [];
toDelete = n_theta_A_no_FS.Newindex ==0;
n_theta_A_no_FS(toDelete,:) = [];
toDelete = n_theta_B_no_FS.Newindex ==0;
n_theta_B_no_FS(toDelete,:) = [];
toDelete = n_beta_A_FS.Newindex ==0;
n_beta_A_FS(toDelete,:) = [];
toDelete = n_beta_B_FS.Newindex ==0;
n_beta_B_FS(toDelete,:) = [];
toDelete = n_beta_A_no_FS.Newindex ==0;
n_beta_A_no_FS(toDelete,:) = [];
toDelete = n_beta_B_no_FS.Newindex ==0;
n_beta_B_no_FS(toDelete,:) = [];


%% 


figure
foo1 = mean(reshape(s_theta_A_FS.Phase,35000,length(range_FS)));
foo2 = mean(abs(reshape(s_theta_A_FS.Spikes,35000,length(range_FS))));
plot(foo1,foo2);
hold on;
errorghost(abs(reshape(s_theta_A_FS.Spikes,35000,length(range_FS))),foo1,'b');

foo1 = mean(reshape(s_theta_A_no_FS.Phase,35000,length(range_FS)));
foo2 = mean(abs(reshape(s_theta_A_no_FS.Spikes,35000,length(range_FS))));
plot(foo1,foo2);
errorghost(abs(reshape(s_theta_A_no_FS.Spikes,35000,length(range_FS))),foo1,'r');
title('Theta dependence of assembly A');

figure
foo1 = mean(reshape(s_theta_B_FS.Phase,35000,length(range_FS)));
foo2 = mean(reshape(s_theta_B_FS.Spikes,35000,length(range_FS)));
plot(foo1,foo2);
hold on;
errorghost(reshape(s_theta_B_FS.Spikes,35000,length(range_FS)),foo1,'b');

foo1 = mean(reshape(s_theta_B_no_FS.Phase,35000,length(range_FS)));
foo2 = mean(reshape(s_theta_B_no_FS.Spikes,35000,length(range_FS)));
plot(foo1,foo2);
errorghost(reshape(s_theta_B_no_FS.Spikes,35000,length(range_FS)),foo1,'r');
title('Theta dependence of assembly B');

% figure
% foo1 = mean(reshape(s_beta_A_FS.Phase,35000,length(range_FS)));
% foo2 = mean(abs(reshape(s_beta_A_FS.Newindex,35000,length(range_FS))));
% plot(foo1,foo2);
% hold on;
% errorghost(abs(reshape(s_beta_A_FS.Newindex,35000,length(range_FS))),foo1,'b');
% 
% foo1 = mean(reshape(s_beta_A_no_FS.Phase,35000,length(range_FS)));
% foo2 = mean(abs(reshape(s_beta_A_no_FS.Newindex,35000,length(range_FS))));
% plot(foo1,foo2);
% errorghost(abs(reshape(s_beta_A_no_FS.Newindex,35000,length(range_FS))),foo1,'r');
% title('Beta dependence of assembly A');

figure
foo1 = movmean(n_beta_A_FS.Phase,10000);
foo2 = movmean(n_beta_A_FS.Newindex,10000);
plot(foo1,foo2);
hold on;
top = foo2+movstd(n_beta_A_FS.Newindex,10000);
bot = foo2-movstd(n_beta_A_FS.Newindex,10000);
plot(foo1,top,'b')
plot(foo1,bot,'b')

foo1 = movmean(n_beta_A_no_FS.Phase,10000);
foo2 = movmean(n_beta_A_no_FS.Newindex,10000);
plot(foo1,foo2);
top = foo2+movstd(n_beta_A_no_FS.Newindex,10000);
bot = foo2-movstd(n_beta_A_no_FS.Newindex,10000);
plot(foo1,top,'r')
plot(foo1,bot,'r')
title('Beta dependence of assembly A');

figure
foo1 = movmean(n_beta_B_FS.Phase,10000);
foo2 = movmean(n_beta_B_FS.Newindex,10000);
plot(foo1,foo2);
hold on;
top = foo2+movstd(n_beta_B_FS.Newindex,10000);
bot = foo2-movstd(n_beta_B_FS.Newindex,10000);
plot(foo1,top,'b')
plot(foo1,bot,'b')
%pgon = polyshape([foo1 fliplr(foo1)],[top fliplr(bot)],'Simplify', false);
%b = plot(pgon,'HandleVisibility','off', 'FaceColor','b');
%b.EdgeAlpha = 0;
%b.FaceAlpha = .15;

foo1 = movmean(n_beta_B_no_FS.Phase,10000);
foo2 = movmean(n_beta_B_no_FS.Newindex,10000);
plot(foo1,foo2);
top = foo2+movstd(n_beta_B_no_FS.Newindex,10000);
bot = foo2-movstd(n_beta_B_no_FS.Newindex,10000);
plot(foo1,top,'r')
plot(foo1,bot,'r')
%pgon = polyshape([foo1 fliplr(foo1)],[top fliplr(bot)],'Simplify', false);
%b = plot(pgon,'HandleVisibility','off', 'FaceColor','r');
%b.EdgeAlpha = 0;
%b.FaceAlpha = .15;
title('Beta dependence of assembly B');

figure
imagesc(new_index_A_FS)
colorbar
figure
imagesc(new_index_A_no_FS)
colorbar
figure
imagesc(new_index_B_FS)
colorbar
figure
imagesc(new_index_B_no_FS)
colorbar

figure
imagesc(spikes_in_A_FS)
figure
imagesc(spikes_in_A_no_FS)
figure
imagesc(spikes_in_B_FS)
figure
imagesc(spikes_in_B_no_FS)

% figure
% plot(s_beta_A_FS.Phase,s_beta_A_FS.Spikes);
% hold on;
% plot(s_beta_A_no_FS.Phase,s_beta_A_no_FS.Spikes);
% 
% figure
% plot(s_beta_B_FS.Phase,s_beta_B_FS.Spikes);
% hold on;
% plot(s_beta_B_no_FS.Phase,s_beta_B_no_FS.Spikes);
