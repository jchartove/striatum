%%  MA666 - Lab 2 (10):  CFC and bicoherence
%   In this tutorial we will first develop our own code to assess cross
%   frequency coupling and bicoherence within a signal.  We will also apply
%   these measures to two case studies of data.  Please note that there 
%   are lots of concepts to consider here.  We'll gloss over most of them,
%   with the ultimate aim of building a CFC measure.  There are many topics to
%   further explore on your own (e.g., filtering, Hilbert transforms, statistical
%   characteristics of measures) if you're interested.

%%  Preliminaries.
%   Text preceded by a '%' indicates a 'comment'.  This text should appear
%   green on the screen.  I will use comments to explain what we're doing 
%   and to ask you questions.  Also, comments are useful in your own code
%   to note what you've done (so it makes sense when you return to the code
%   in the future).  It's a good habit to *always* comment your code.  I'll
%   try to set a good example, but won't always . . . 

%%  Part 1:  A Hilbert transform example.
%   An important component of the CFC measure discussed in class is the
%   determination of the phase or amplitude envelope of a signal.  To do
%   so, we'll compute the Hilbert transform.
%
%   Let's start by considering a single sinusoid.

dt = 0.001;             %Sampling interval [s].
N = 1000;               %Number of samples.
t = (0:N-1)*dt;         %Time axis [s].
f1 = 2.0;               %Freq of sinusoid [Hz].
phi0 = 0.0;             %Initial phase of sinusoid.
d = sin(2.0*pi*t*f1 + phi0);

%   And compute the analytic signal,

dA = hilbert(d);

%   Note that the call to 'hilbert' in MATLAB returns the analytic signal!
%   We can plot the original signal, and the real and imaginary parts of
%   the analytic signal.

figure(10);  clf()
subplot(4,1,1)
plot(d);  ylabel('Data');
subplot(4,1,2)
plot(real(dA))
hold on
plot(imag(dA), 'r')
hold off
ylabel('Re (blue), Im (red)')
axis tight

%   Notice that the imaginary part of the analytic signal is the real part
%   of the analytic signal shifted by 90 degress (i.e., by pi/2).

%   Finally, we can compute the phase of the analytic signal,

phi = angle(dA);

%   and the amplitude envelope,

amp = abs(dA);

%   and let's plot both results.

subplot(4,1,3)
plot(phi)
ylabel('Angle')
axis tight

subplot(4,1,4)
plot(amp)
ylabel('Amplitude')
axis tight

%   Notice that the angle increases linearly from -pi to pi, and at pi
%   returns to -pi.  Also notice that, in this case, the amplitude of the
%   signal is one for all time.  Make sense?

%%  Part 2:  The analytic signal has no power at negative frequencies.
%   In class, we discussed that the analytic signal has no power at negative frequencies
%   and we showed this for a very simple example.  Let's verify this result
%   here in a numerical simulation.

%   First, define a signal.  We'll make it noise.
dt = 0.001;  N=1000; t=(0:N-1)*dt;  T=N*dt;
d = randn(N,1);

%   Then, compute the analytic signal.
dA = hilbert(d);

%   Now, compute the power spectra of the orginal signal and the
%   analytic signal.  Do this for *all* frequencies, both positive and
%   negative.

pow_d  = 2*dt^2/T * (fft(d) .*conj(fft(d)));
pow_dA = 2*dt^2/T * (fft(dA).*conj(fft(dA)));

%   Reorganize the power spectra to start at negative frequencies and
%   increases through zero and then positive frequencies,

pow_d  = fftshift(pow_d);
pow_dA = fftshift(pow_dA);

%   Define the frequency axis, and include both postive and negative freq.
df = 1/T;
fNQ = 1/dt/2;
faxis = (-fNQ:df:fNQ-df);

figure(10);  clf()
subplot(2,1,1)
plot(t,d);  ylabel('Data');  xlabel('Time [s]')

subplot(2,1,2)
plot(faxis, pow_d);  ylabel('Power of Signal');  xlabel('Freq [Hz]')
hold on
plot(faxis, pow_dA/4, 'r');  title('Power of Analytic Signal (red)')
hold off

%   Notice that i) For the analytic signal, the negative frequencies have zero
%   power, and ii) We scale the power of the analytic signal by 1/4 to match the
%   power of the original signal.  The reason:  the analytic signal doubles
%   the amplitude at positive frequencies, and therefore increase the power by
%   a factor of 4 (because power is amplitude squared).
%
%   We conclude that the analytic signal is the original signal with no power
%   at negative frequencies.

%%  Part 3:  CFC.
%   Now, let's compute the cross-frequency coupling (CFC) of a signal.  To
%   start, let's load data with true cross-frequency coupling.  You'll need
%   to download the file,
%
%       data_1.mat
%
%   from the Blackboard.  Then, load these data into MATLAB,

clear                               %First, clear the workspace,
load('data_1.mat')                  %... then load the data.

%   Plot the data, and check through visual inspection whether CFC occurs,

N = length(d);  Fs = 1/dt;          %Define useful parameters.
t = (0:N-1)*dt;                     %Define a time axis.
plot(t, d)
xlabel('Time [s]')

%IN LAB Q:  Can you see any CFC?

%   To assess the rhythmic activity, compute the  power spectrum,

periodogram(d, [], [], 1/dt)

%   From this figure, we observe two peaks.  i) A sharp low frequency peak near 6 Hz,
%   and ii) A broad high frequency peak from 60-140 Hz.  These two regions
%   will define our low and high frequency bands.  Notice that the high
%   frequency activity has much lower power than the low frequency activity.
%
%   Now, we'll develop a measure to assess the cross frequency coupling.  To do
%   so, we'll implement three steps discussed in class.
%
%   Step 1:  Filter the data into low and high frequency bands.  To filter
%   the data, we'll use a second order Butterworth filter.  (Note that
%   there are many possible choices here!)

deg=2;                              %Sets the filter order.
Wn = [5*2/Fs,7*2/Fs];               %Define the low frequency window of interest.
[B,A] = butter(deg,Wn,'bandpass');  %Apply the filter to isolate the band.
dlo = filtfilt(B,A,d);

Wn = [60*2/Fs,140*2/Fs];            %Define the high frequency window of interest.
[B,A] = butter(deg,Wn,'bandpass');  %Apply the filter to isolate the band.
dhi = filtfilt(B,A,d);

%   How well does the filter work?  Let's compare the original signal with
%   the fitlered data.

plot(t,d)
hold on
plot(t,dlo, 'r')
plot(t,dhi, 'g')
hold off
xlabel('Time [s]');  title('Data (blue), Low-freq (red), High-freq (green)')

%   We might notice that funny effects happen at the edges of the data.  So,
%   let's ignore these edge effects due to the filtering and concentrate on
%   the center of the data.

dlo = dlo(N/4:end-N/4-1);
dhi = dhi(N/4:end-N/4-1);
taxis = t(N/4:end-N/4-1);

%   Step 2:  Compute the phase of the low frequency signal and the
%   amplitude envelope of the high frequency signal.  To do so, we use the
%   Hilbert transform.

phi = angle(hilbert(dlo));    %Phase of low frequency signal.
amp = abs(hilbert(dhi));      %Amplitude envelope of high frequency signal.

%   Let's check that these results look okay,

figure(10);  clf()
subplot(2,1,1)
plot(taxis, dlo)
hold on
plot(taxis, phi, 'g')
hold off
axis tight
xlabel('Time [s]');  title('Low-freq and phase');

subplot(2,1,2)
plot(taxis, dhi)
hold on
plot(taxis, amp, 'r')
hold off
axis tight
xlabel('Time [s]');  title('High-freq and amplitude envelope')

%   Step 3:  Determine if the phase and envelope are related.  To do so,
%   we'll divide the phases into bins, find the times where the low
%   frequency phase lies in each bin, and then compute the average
%   amplitude envelope of the high frequency signal at those times.

p_bins = (-pi:0.2:pi);              %Define the phase bins.

a_mean = zeros(length(p_bins)-1,1); %Vector to hold avg amp env results.
p_mean = zeros(length(p_bins)-1,1); %Vector to hold center of phase bins.
  
for k=1:length(p_bins)-1
    pL = p_bins(k);                         %Phase lower limit for this bin.
    pR = p_bins(k+1);                       %Phase upper limit for this bin.
    indices = find(phi >= pL & phi < pR);   %Find phase values in this range.
    a_mean(k) = mean(amp(indices));         %Compute mean amplitude at these phases.
    p_mean(k) = mean([pL, pR]);             %Label the phase bin with the center phase.
end
h = max(a_mean)-min(a_mean);                %The difference between the max and min modulation.

%   Plot the mean envelope versus phase.
figure(11)
plot(p_mean, a_mean, 'k', 'LineWidth', 1);
axis tight
xlabel('Low frequency phase');  ylabel('High frequency envelope height difference');
title(['Metric h=' num2str(h)])

%   In this figure, the relationship between the phase and amplitude
%   becomes apparent.  So, we conclude CFC exists in these data.

%%  Part 4:  Bicoherence.
%   The second measure of coupling we'll consider is the
%   bicoherence.  To have bicoherence in a signal, we require a
%   constant phase relationship (over trials) between oscillations
%   at three frequencies, f1, f2, and their sum f1+f2.  To compute
%   the bicoherence, download the function from Blackboard,
%
%   bicoherence.m
%
%   In class, we'll look at this code and see if it makes sense.
%
%   To explore how this measure works, we'll create some simulated
%   data and compute the bicoherence.  Let's start with pure noise,

clear

K=100; dt=0.001; N=1000;
x = zeros(K,N);
 
for k=1:K
    x(k,:) = randn(1,N);
end

%   Notice that we create multiple trials of the data.  To compute
%   the bicoherence,

[b, faxis] = bicoherence(x, 1/dt, 50);

%   Notice that the result is two dimensional.  To display the
%   results,

imagesc(faxis, faxis, b, [0,1])
colorbar
axis xy
xlabel('Frequency [Hz]');  ylabel('Frequency [Hz]');  title('Bicoherence');

%IN LAB Q: What do you find for the bicoherence?  Does it match your
%expectation?