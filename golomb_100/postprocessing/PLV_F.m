%stim = decimate(model.fixed_variables.dend_iPeriodicPulsesBen_s2(:,1),10);
%foi = model.parameters.dend_iPeriodicPulsesBen_PPfreq;
%[plv, plv_s] = PLV_F(soma_V, stim, foi, time);

function [ plv, plv_s] = PLV_F(signal, stim, foi, time)
% Compute the Phase Locking Value between two signals across trials, according to Lachaux,
% Rodriguez, Martinerie, and Varela (1999). The PLV value ranges from 0, indicating random
% phase differences, to 1 indicating a fixed phase difference.
% phase_sig1 and phase_sig2 should be the phase values of the signals in radians, arranged as
% Samples x Trials. These can bed
% computed using the Wavelet or Hilbert transform, for example:
% phase_sig = angle(hilbert(BPS));
% Where BPS is the signal after band-pass filtering around the frequency range of interest.
%
% JC edit 4/1/20: building in the phase computation. signal should be soma_soma_somaSomaiSYN_s.
% foi means frequency of interest
%
% Written by Edden Gerber 2012

dt = 0.01; %hardcoding for personal convenience

%preprocessing
n = size(signal); %check dimensionality
if n(1) > 1
    mean_signal = nanmean(signal, 2);
else
    mean_signal = signal;
end
low_time = 500;
time_index = time >= low_time;
mean_signal_detrended = detrend(mean_signal(time_index, :));

%[s,f,t] = spectrogram(x,window,noverlap,f,fs) returns the spectrogram at the cyclical frequencies specified in f.
[s,w,t1] = spectrogram(mean_signal_detrended,1200,1100,[0:200],100/dt,'yaxis'); %sampling frequency is decimated
BPS = s(w == foi, :);
phase_sig = angle(hilbert(BPS));
t1_adjusted = t1*1000 + 500;
time2 = logical(sum(time==t1_adjusted,2));

%convert square wave to sinusoid
stim_sin = wavelet_spectrogram(stim, 10000, foi, foi, 1); %this is expensive but the easiest way to do this
stim_phase = angle(hilbert(stim_sin));
stim_phase_short = stim_phase(time_index,:);

%for spikes
T_short = (size(mean_signal_detrended,1)-1);
numcells = size(signal,2);
V_short = signal(time_index,:);
spike_indicator = zeros(numcells,T_short);
spikes = zeros(1,numcells);

for t = 1:T_short
    spike_indicator(:,t) = (V_short(t,:)<0) & (V_short(t+1,:) >= 0);
    s = (V_short(t,:)<0) & (V_short(t+1,:) >= 0);
    spikes = spikes + s;
end

%[s,f,t] = spectrogram(x,window,noverlap,f,fs) returns the spectrogram at the cyclical frequencies specified in f.
%[s,w,t] = spectrogram(sum(spike_indicator),1200,1100,[0:100],100/dt,'yaxis'); %sampling frequency is decimated
%BPS_spikes = s(w == foi, :);

%should be phase of input at spike time (sum of stimulus phase evaluated at spike indicator)
stim_phase_short=stim_phase_short(1:size(spike_indicator,2));
phase_spikes = logical(spike_indicator).*stim_phase_short';

% compute PLV
[~, Ntrials] = size(phase_sig);
e = exp(1i*(phase_sig - stim_phase(time2)'));
plv = abs(sum(e,2)) / Ntrials;

% compute PLV_s
plv_s = abs(nanmean(exp(sqrt(-1).*phase_spikes(find(phase_spikes)))));
end
