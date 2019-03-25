%%  MA666 - Lab 3 (11):  Spike-field coherence. 
%   In this tutorial we will develop tools to compute the spike-field
%   coherence.  We'll first see how to implement this code ourselves, then
%   we'll interface with the MATLAB software package Chronux. 

%%  Preliminaries.
%   Text preceded by a '%' indicates a 'comment'.  This text should appear
%   green on the screen.  I will use comments to explain what we're doing 
%   and to ask you questions.  Also, comments are useful in your own code
%   to note what you've done (so it makes sense when you return to the code
%   in the future).  It's a good habit to *always* comment your code.  I'll
%   try to set a good example, but won't always . . . 

%%  Part 1:  Power spectrum of a spike train.
%   To compute the power spectrum of a spike train, we first need to calculate
%   the Fourier transform of the spike train.  We'll start by considering
%   the case of random (Poisson) spiking, and compute the Fourier transform
%   (FT) both 'by-hand' and using built-in MATLAB routines.  We already can
%   guess what this spectrum should look like, right?

%   Let's begin by clearing the workspace.
clear

%   First, we'll create simulated data.
dt = 0.001;             %The sampling interval, 1 ms in this case.
N = 1000;               %The number of time steps.
t = (0:N-1)*dt;         %Define a time axis.
T = N*dt;               %Define the total time of the simulation.

%   Now, let's set the firing rate, lambda, which indicates the expected
%   number of spikes per second,
lambda = 10;            

%   To generate Poisson spiking, we'll generate N Poisson random numbers.
%   Note that the probability of a spike in any (tiny) bin is the firing
%   rate lambda multiplied by the size of the bin (dt),

dn = poissrnd(lambda*dt,1,N);    %The simulated spike train.

%   Now, let's subtract the mean rate from the spike train.
dn = dn - mean(dn);

%   And let's plot it.
subplot(3,1,1)
plot(t, dn)
xlabel('Time [s]');  ylabel('Spikes')

%% IN LAB Q:  How many spikes to you observe per second?  Is it consistent
%with lambda?

%   Before we compute the power spectrum, let's first compute the
%   auto-covariance of the spike train.  To do so, we'll use the built-in
%   MATLAB function xcorr.  Can you guess what the auto-covariance will look
%   like?

[acf, lag] = xcorr(dn, 'biased');       %Compute the AC.
acf = acf(N/2:end-N/2);                 %Only keep the middle half of AC.
lag = lag(N/2:end-N/2);                 %... and of the lags.

%   Let's also plot the AC.  Does it make sense?
subplot(3,1,2)
plot(lag, acf)

%   Finally, we'll compute the FT of the spike train.  We'll first do so
%   by-hand.  To start, define the frequency axis.  Here we arrange the
%   frequencies following the MATLAB convention.
fj = [(0:N/2-1)*(1/T), (-N/2:-1)*(1/T)];

%   We can also define a taper - here we use the 'default' taper.
h = ones(1,N);

%   Then compute the FT 'by-hand'.
X0 = zeros(1,length(fj));                           %An empty vec to hold results.
for j=1:length(fj)                                  %For each freq,
    X0(j) = sum(h .* dn .* exp(-2*pi*1i*fj(j)*t));  %... data * taper * sinusoids.
end

%   The power spectrum is the FT multiplied by its complex conjugate, and
%   scaled.  Let's plot the results.
subplot(3,1,3)
S = dt^2/T*(X0.*conj(X0));                       %Compute the power.
plot(fftshift(fj),fftshift(S))                   %... and plot it, with increasing freq. 

%   Above we computed the FT by-hand.  We don't need to do this.  Instead,
%   we can compute the power spectrum using built-in MATLAB routines.  Do
%   so, and plot the results, in the next few lines.
pow = dt^2/T*(abs(fft(dn)).^2);
hold on
plot(fftshift(fj), fftshift(pow), 'og')
hold off

%TECHNICAL NOTE:  
%In addition, we can also plot the *expected* power spectrum for the Poisson
%process.  We said in class that, for random spiking, we expect a
%flat power spectrum (i.e., all frequencies to be present).  We
%stated that the power of this flat spectrum would equal "lambda*dt".
%In fact, we need to plot "dt" times this expected quantity.  In
%this case, we've used the same normalization of the power spectrum
%introduced in MA665, namely the term "dt^2/T".  By doing so, we need to
%include an extra factor of "dt" when plotting the expected spectrum.

hold on
plot(fj,dt*(dt*lambda)*ones(length(fj),1), 'r')  %Also, plot the expected spectrum.
hold off

%%  Part 2:  Power spectrum for multiple trials of Poisson spiking data.
%   Above, we computed the power spectrum of a single spike train.  It's
%   often the case that multiple instances of a spike train will be
%   observed (over multiple trials, say).  Let's consider multi-trial
%   observations of a spike train, and compute the average spectrum over
%   the trials.

%   Let's begin by clearing the workspace.
clear

%   Then, define the number of trials, the number of points in each trial, the time
%   step, and the total time of each trial.
K=20; N = 1000; dt = 0.001; T=N*dt;

%   Then define the frequencies for the power spectrum.
fj = [(0:N/2-1)*(1/T), (-N/2:-1)*(1/T)];

%   We'll again consider Poission spikes with fixed rate lambda.
lambda = 10;

%   ... and generate multiple trials.  For each trial, generate the Poisson
%   spike data, compute the power spectrum, and compute the AC.  Save the
%   results of these last two computations for each trial.
S = zeros(K,length(fj));
acf = zeros(K, N);
for k=1:K

    dn = poissrnd(lambda*dt,1,N);
    dn = dn - mean(dn);

    X0 = fft(dn);
    S0 = dt^2/T*(abs(fft(dn)).^2);
    
    S(k,:) = S0;

    [acf0, lag] = xcorr(dn, 'biased');       %Compute the AC.
    acf(k,:) = acf0(N/2:end-N/2);           %Only keep the middle half of AC.
end
lag = lag(N/2:end-N/2);

%   Finally, plot the results of the AC and the power spectrum, both
%   averaged over all K trials.

subplot(2,1,1)
plot(lag, mean(acf,1)); xlabel('Lags [ms]'); ylabel('AC')

subplot(2,1,2)
S = mean(S,1);
plot(fftshift(fj),fftshift(S))                 %... and plot it, with increasing freq. 
hold on
plot(fj,dt*(dt*lambda)*ones(length(fj),1), 'r')
hold off
ylim([0 0.0001]); xlabel('Freq [Hz]'); ylabel('Power')

%%  Part 3:  Power spectrum for multiple trials of refractory spiking data.
%   So far, we only considered examples of Poisson random data.  We did so
%   because we knew what power spectrum to expect (a constant one).  Let's
%   now consider a second example of spike train data:  spiking with a
%   refractory period.  We discussed in lecture what the spectrum should
%   be.  Let's see if that expectation holds.

%   The code below is the same as the code in Part 2.  What's changed is
%   the data we'll analyze.

%   Let's begin by clearing the workspace.
clear

%  Initial set up of simulation parameters.
K=50; N=1000; dt=0.001; T=N*dt;

%  Define the frequencies.
fj = [(0:N/2-1)*(1/T), (-N/2:-1)*(1/T)];

%  Define quantities that change over trials.
lambda=zeros(K,1);
S = zeros(K,length(fj));
acf = zeros(K, N);
for k=1:K

    %Generate Poisson spiking with a refractory period.  This part is a bit
    %tricky.  In the above example, we fixed lambda - the probability of a 
    %spike.  Here, we allow lambda to vary in time.  Specifically, after a
    %spike occurs, we make the probability of a spike in the next 5 steps
    %very small;  this acts as the refractory period.
    c0=100;  c1=100;  b=-0.5;
    dn = zeros(1,N);
    lambda0 = zeros(1,N);
    for t=6:N
        lambda0(t) = c0 + c1*exp(b*dn(t-1)+b*dn(t-2)+b*dn(t-3)+b*dn(t-4)+b*dn(t-5));
        dn(t) = poissrnd(lambda0(t)*dt); 
    end
    lambda(k) = mean(lambda0);
    
    %The rest of the code is the same as Part 2.
    
    dn = dn - mean(dn);
    
    X0 = fft(dn);
    S0 = dt^2/T*(abs(fft(dn)).^2);
    
    S(k,:) = S0;

    [acf0, lag] = xcorr(dn, 'biased');       %Compute the AC.
    acf(k,:) = acf0(N/2:end-N/2);
end
lag = lag(N/2:end-N/2);

%   And plot the results.
subplot(2,1,1)
plot(lag, mean(acf,1));
xlabel('Lag [ms]');  ylabel('AC')

subplot(2,1,2)
S = mean(S,1);
plot(fftshift(fj),fftshift(S))
hold on
plot(fj,dt*(dt*mean(lambda))*ones(length(fj),1), 'r')
hold off
ylim([0 0.0003]);  xlabel('Freq [Hz]'); ylabel('Power')

%IN LAB Q:  Are the AC and the power spectrum as you expect?

%%  Part 4:  Spike-spike coherence. 
%   We've now computed the FT of spike trains, and used the FT to compute
%   the power spectrum.  We can also use the FT to compute the coherence
%   between two spike trains.  Remember that coherence detects a constant
%   phase relationship between two siganls over trials at some frequency.
%   The multi-trial structure is critical (just as it was for the field
%   data we considered earlier in the course).
%
%   We'll start by considering two sets of Poisson spike trains.  We
%   already have an expectation for the coherence (right?).  Let's now
%   compute it.

%   Let's begin by clearing the workspace.
clear

%   Set the simulation parameters, using multiple trials.
K=20; N = 1000; dt = 0.001; T=N*dt;
t = (0:N-1)*dt;

%   Define the frequencies.
fj = [(0:N/2-1)*(1/T), (-N/2:-1)*(1/T)];

%   We'll consider random spikes with rate lambda.
lambda = 10;

%   Define output variables for the Fourier transform of the two spike trains.
X1 = zeros(K,length(fj));
X2 = zeros(K,length(fj));
for k=1:K

    %Generate the random data.
    dn1 = poissrnd(lambda*dt,1,N);
    dn2 = poissrnd(lambda*dt,1,N);

    %Compute mean rate of each.
    lambda1 = mean(dn1);
    lambda2 = mean(dn2);

    %Subtract the mean rates.
    dn1 = dn1-mean(dn1);
    dn2 = dn2-mean(dn2);
    
    %Compute the FTs.
    X01 = fft(dn1);
    X02 = fft(dn2);
    
    %Save the results.
    X1(k,:) = X01;
    X2(k,:) = X02;
end

%   Now, compute the power spectra and cross spectrum.
S11 = zeros(1,length(fj));
S12 = zeros(1,length(fj));
S22 = zeros(1,length(fj));
for k=1:K
    S11 = S11 + dt^2/T*(X1(k,:).*conj(X1(k,:)))/K;
    S12 = S12 + dt^2/T*(X1(k,:).*conj(X2(k,:)))/K;
    S22 = S22 + dt^2/T*(X2(k,:).*conj(X2(k,:)))/K;
end

%   And then the coherence.
cohr = S12.*conj(S12) ./S11 ./S22;

%   Let's plot everything to have a look.
subplot(3,1,1)
plot(t, dn1)
xlabel('Time [s]'); ylabel('Neuron 1')

subplot(3,1,2)
plot(t, dn2)
xlabel('Time [s]'); ylabel('Neuron 2')

subplot(3,1,3)
plot(fftshift(fj),fftshift(cohr))
ylim([0 1])
xlabel('Freq [Hz]');  ylabel('Coherence')

%IN LAB Q: Are the results what you expected?

%%  Part 5:  Spike-field coherence.
%   In the previous section we considered the coherence between two spike
%   trains.  Now let's consider the coherence between a spike train, and a
%   field (for example, an LFP).  We still require a multi-trial structure.
%   The code that follows is very similar to the code in Part 4;  the
%   primary change is to replace one of the spike trains with a field.

%   Let's begin by clearing the workspace.
clear

%   Set the simulation parameters, using multiple trials.
K=20; N = 1000; dt = 0.001; T=N*dt;
t = (0:N-1)*dt;

%   Define the frequencies.
fj = [(0:N/2-1)*(1/T), (-N/2:-1)*(1/T)];

%   We'll include random spikes with rate lambda0.
lambda = 10;

%   We'll also generate an LFP with frequency f1 Hz.
f1 = 10;

%   Define output variables for the spectra of the two signals.
X1 = zeros(K,length(fj));
X2 = zeros(K,length(fj));
for k=1:K

    %Define the LFP signal as a sinusoid + noise.
    s  = sin(2.0*pi*t * f1 + 1*pi*randn)+0.1*randn(1,N);
    
    %The spikes occur at random times.
    dn = poissrnd(lambda*dt,1,N);

    %Subtract the mean rates.
    dn = dn-mean(dn);
    
    %Compute the FTs.
    X01 = fft(dn);
    X02 = fft(s);
    
    %Save the results.
    X1(k,:) = X01;
    X2(k,:) = X02;
end

%   Now, compute the power spectra and cross spectrum.
S11 = zeros(1,length(fj));
S12 = zeros(1,length(fj));
S22 = zeros(1,length(fj));
for k=1:K
    S11 = S11 + dt^2/T*(X1(k,:).*conj(X1(k,:)))/K;
    S12 = S12 + dt^2/T*(X1(k,:).*conj(X2(k,:)))/K;
    S22 = S22 + dt^2/T*(X2(k,:).*conj(X2(k,:)))/K;
end

%   And then the coherence.
cohr = S12.*conj(S12) ./S11 ./S22;

%   Let's plot everything to have a look.
subplot(3,1,1)
plot(t, dn)
xlabel('Time [s]'); ylabel('Neuron 1')

subplot(3,1,2)
plot(t, s);
xlabel('Time [s]'); ylabel('Field')

subplot(3,1,3)
plot(fftshift(fj),fftshift(cohr))
ylim([0 1]);  xlim([-50 50])
xlabel('Freq [Hz]');  ylabel('Coherence')

%IN LAB Q: Are the results what you expected?

%%  Part 6:  Chronux.
%   Above we computed the coherence (between two spike trains, or between a
%   spike train and a field) "by hand". There's a nice software package
%   for doing this, which you can download from here,
%
%   http://www.chronux.org
%
%   Once you've installed this software, it's easy to compute the coherence
%   between two sets of spike trains.  This software implements what's
%   called the "multi-taper method".  The idea is that - instead of
%   applying a single taper to the data - we apply multiple tapers, and
%   then average the results across tapers.  This averaging reduces the
%   noise in our spectral estimates.  But, there's a cost:  reduced
%   frequency resolution.  There are many details to explore here;  if you're
%   interested, check out this week's chapter, and Chapter 5, section 3.8.
%   An example of generating artificial data (the same as in Part 5) and
%   applying the Chronux software is provided below.

%   Let's begin by clearing the workspace.
clear

%   Set the simulation parameters, using multiple trials.
K=50; N = 1000; dt = 0.001; T=N*dt;
t = (0:N-1)*dt;

%   We'll include random spikes with rate lambda0.
lambda = 10;

%   We'll also generate an LFP with frequency f1 Hz.
f1 = 10;

%   Define multi-trial spike train and field data.
s = zeros(K,N);
n = zeros(K,N);
for k=1:K

    %Define the LFP signal as a sinusoid + noise.
    s(k,:)  = sin(2.0*pi*t * f1 + 1*pi*randn)+0.1*randn(1,N);
    
    %The spikes occur at random times.
    dn(k,:) = poissrnd(lambda*dt,1,N);

end

%%%%%%%%  Here's the Chronux part!  %%%%%%%%

%   Set up parameters for Chronux.
params.Fs = 1/dt;       %Set the sampling frequency.
params.tapers = [5 9];  %Set the time-bandwidth product TW and # tapers.
params.trialave = 1;    %Average the results over trials. 

%   Compute the coherence.
[C,phi,S12,S1,S2,f]=coherencycpb(s',dn',params);

%   Plot the results.
plot(f,C)
xlabel('Freq [Hz]');  ylabel('Coherence')
ylim([0 1])

