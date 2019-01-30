%  Due: Nov 25, 2014
%
%  x.  Please create a PDF with your answers.
%  x.  Please email you answers to:  mak@bu.edu
%  x.  In your email, please use the Subject:  MA666 Lab 3
%  
%  x.  Please answer at least **three** "Challenge Problems".
%  Your solutions should include figures (with axes labeled
%  and captions) and explanatory text.  Please include any
%  MATLAB code you write;  I'll be looking for comments in your
%  code!  The comments should be detailed enough that you could share your
%  code with any of your classmates / colleagues, and he/she could easily
%  read and use your code.
%
%  x.  You may work together.

%% Challenge Problem #1.
%Analyze your own data set of spiking neurons, or of spiking neurons and
%fields, using the tools developed in this lab.  Deduce a scientific
%conclusion.

%% Challenge Problem #2.
%In Part 5 of Lab, we considered spike-field coherence between
%a sinusoidal field and random spike train, and found little spike-field
%coherence, as expected.  In this Challenge, please generate a simulated data
%set with *strong* spike-field coherence.  You might start with the rhythmic field
%('s') in Part 5, and adjust the spiking to occur at a particular phase of
%the field.  Compute the spike-field coherence, and show that it is near 1
%at some frequency, and near zero at other frequencies.

%% Challenge Problem #3.
%Consider the data set 'case_study_1.mat' on Blackboard for this week's lab.
%The data consists of three components,
%
%    s = An LFP field.
%    n1 = A spike train from Neuron #1.
%    n2 = A spike train from Neuron #2.
%
%Each component is 100 trials by 1000 indices, and the sampling interval
%between indices (i.e., dt) is 1 ms.  Analyze these data using the tools
%developed in this lab.  Specifically, consider:

plot(t,mean(n1),'LineWidth',2)
hold on
plot(t,mean(n2),'r')
xlabel('Time (s)')
ylabel('Amplitude')
plot(t,mean(s),'g', 'LineWidth',2)
title('Mean data, Case 1 (n1 in blue, n2 in red, s in green)')
%
%(Q.i) The field power spectrum, averaged across trials.  Does the field
%have a dominant rhythm?
dt = 0.001;             %The sampling interval, 1 ms in this case.
N = 1000;               %The number of time steps.
t = (0:N-1)*dt;         %Define a time axis.
T = N*dt;               %Define the total time of the simulation.
K = 100;

fj = [(0:N/2-1)*(1/T), (-N/2:-1)*(1/T)];
h = ones(1,N);

powF = zeros(K,length(fj));
for k=1:K
    X = fft(s(k,:));
    powF(k,:) = powF(k,:) + 2*dt^2/T* X.*conj(X);
end

powF = sum(powF)/K;
plot(fj, 10*log10(powF));
xlabel('Freq [Hz]');  ylabel('Power [dB]');
xlim([0 50])
title('Field power spectrum, Case 1')

%
%(Q.ii) The spike power spectra, averaged across trials.  Do the spike
%trains from Neuron #1 and Neuron #2 have a dominant rhythm?
%

pow1 = zeros(K,length(fj));
pow2 = zeros(K,length(fj));
for k=1:K

    dn1 = n1(k,:);
    dn1 = dn1 - mean(dn1);
    
    dn2 = n2(k,:);
    dn2 = dn2 - mean(dn2);

    X01 = fft(dn1);
    S01 = dt^2/T*(abs(fft(dn1)).^2);
    
    X02 = fft(dn2);
    S02 = dt^2/T*(abs(fft(dn2)).^2);
    
    pow1(k,:) = S01;
    pow2(k,:) = S02;
end

pow1 = mean(pow1,1);
pow2 = mean(pow2,1);
plot(fftshift(fj),fftshift(pow1), 'LineWidth', 2);
hold on
plot(fftshift(fj),fftshift(pow2),  'r');
xlabel('Freq [Hz]');  ylabel('Power');
xlim([0 500])
title('Spike power spectrum, Case 1 (n1 in blue, n2 in red)')
xlim([0 50])
title('Spike power spectrum closeup, Case 1 (n1 in blue, n2 in red)')

%(Q.iii) The spike-spike coherence. Are the spike trains coherent?
%
X1 = zeros(K,length(fj));
X2 = zeros(K,length(fj));
for k=1:K
    dn1 = n1(k,:);
    dn2 = n2(k,:);
    
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
plot(fftshift(fj),fftshift(cohr))
ylim([0 1])
xlabel('Freq [Hz]');  ylabel('Coherence')
title('Spike-spike coherence, Case 1')

%(Q.iv) The spike-field coherence. Is each spike train coherent with the field?

%   Now, compute the power spectra and cross spectrum.
X1 = zeros(K,length(fj));
X2 = zeros(K,length(fj));
for k=1:K
    dn1 = n1(k,:);
    dn2 = s(k,:);
    
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

plot(fftshift(fj),fftshift(cohr))
ylim([0 1]);  xlim([-50 50])
xlabel('Freq [Hz]');  ylabel('Coherence')


%% Challenge Problem #4.
%Consider the data set 'case_study_2.mat' on Blackboard for this week's lab.
%The data consists of three components,
%
%    s = An LFP field.
%    n1 = A spike train from Neuron #1.
%    n2 = A spike train from Neuron #2.
%
%Each component is 100 trials by 1000 indices, and the sampling interval
%between indices (i.e., dt) is 1ms.  Analyze these data using the tools
%developed in this lab.  Specifically, consider:
%
%(Q.i) The field power spectrum, averaged across trials.  Does the field
%have a dominant rhythm?
%
%(Q.ii) The spike power spectra, averaged across trials.  Do the spike
%trains from Neuron #1 and Neuron #2 have a dominant rhythm?
%
%(Q.iii) The spike-spike coherence.  Are the spike trains coherent?
%
%(Q.iv) The spike-field coherence. Is each spike train coherent with the
%field?