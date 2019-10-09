%  Due: Nov 6, 2014
%
%  x.  Please create a PDF with your answers.
%  x.  Please email you answers to:  mak@bu.edu
%  x.  In your email, please use the Subject:  MA666 Lab 1
%  
%  x.  Please answer at least the first **three** "Challenge Problems".
%  Your solutions should include figures (with axes labeled
%  and captions) and explanatory text.  Please include any
%  MATLAB code you write;  I'll be looking for comments in your
%  code!  The comments should be detailed enough that you could share your
%  code with any of your classmates / colleagues, and he/she could easily
%  read and use your code.
%
%  x.  You may work together.

%% Challenge Problem #1.
%Update the coherence measure to include multiple trials. Consider the
%following simulated data,
%
   K=100; dt=0.001; N=1000; t=(0:N-1)*dt;
   x = zeros(K,N);
   y = zeros(K,N);

   for k=1:K
     x(k,:) = randn(1,length(t));
     y(k,:) = randn(1,length(t));
   end
%
%Notice that there are 100 trials, each of length 1 s (here dt=1ms).
%Also notice that the data is random in each trial.  So we expect the
%ensemble averaged coherence to be zero.
%
%(Q.i) Update the code developed in Lab to compute the ensemble
%averaged coherence for these pure-noise data.  Plot the resulting
%coherence and compare it to your expectation for what the coherence
%*should* be for these simulated data.
%
%(Q.ii)  Repeat your coherence analysis for the following simulated data
%consisting of two sinusoids with constant phase shift and added noise,
%
    K=100; dt=0.001; N=1000; t=(0:N-1)*dt;
    x = zeros(K,N);  f1 = 10;
    y = zeros(K,N);  f2 = 10;
 
    for k=1:K
      x(k,:) = sin(2.0*pi*t*f1 + 0.3*pi) + 0.1*randn(1,length(t));
      y(k,:) = cos(2.0*pi*t*f1) + 0.1*randn(1,length(t));
    end
%
%Be sure to plot the resulting coherence, and explain your result;  namely,
%does the result make sense with your expectations?

%% Challenge Problem #2.  Case Study 1.
%Download the data 'Case_Study_1.mat' from Blackboard.  You'll find the
%file contains three variables,
%
%  x = observation of time series x [ntrials, time].
%  y = observation of time series y [ntrials, time].
%  dt = sampling interval in [s].
%
%Both x and y are observed simultaneously for 100 trials, each trial lasting 1s.
%The sampling interval dt = 0.001s.  Use the tools developed in this lab
%to analyze these data.  Please compute:

plot(mean(x))
hold on
plot(mean(y),'r')
xlabel('Time (ms)')
ylabel('Amplitude')
title('Average raw data of X (blue) and Y (red), case 1')
%
%(Q.i) The power spectra of x and y.  For these data, compute the power spectra
% for each trial, and average the results over all trials.

N = length(x(1,:));
K = length(x(:,1));
t = (0:N-1)*dt;
T = N*dt; 
fj = [(0:N/2-1)*(1/T), (-N/2:-1)*(1/T)];

powX = zeros(1,N);
powY = zeros(1,N);

for k = 1:K
    X = fft(x(k,:));
    powX = powX + 2*dt^2/T* X.*conj(X);
    Y = fft(y(k,:));
    powY = powY + 2*dt^2/T* Y.*conj(Y);
end

powX = powX/K;
powY = powY/K;
plot(fj, 10*log10(powX), 'LineWidth', 2);
hold on
plot(fj, 10*log10(powY), 'r');
xlabel('Freq [Hz]');  ylabel('Power [dB]');

%(Q.ii) The cross correlation between x and y.  Again, for these data, compute
%the cross correlation for each trial, and average the results over all trials.
%

rxy_sum = zeros(1,2*N-1);
lag_sum = zeros(1,2*N-1);
for k = 1:K
    [rxy, lag] = xcorr(x(k,:),y(k,:),'biased');
    rhoxy = rxy/(var(x(k,:))*var(y(k,:)));
    rxy_sum = rxy_sum + rhoxy;
    lag_sum = lag_sum + lag;
end
rxy_sum = rxy_sum/K;
lag_sum = lag_sum/K;
plot(lag_sum, rxy_sum)
xlabel('Lags [ms]');  ylabel('Cross Correlation');

%(Q.iii) The coherence between x and y.

df = 1/T;
fNQ = 1/dt/2;
faxis = (-fNQ:df:fNQ-df);

sum_Sxx = zeros(1,N);
sum_Syy = zeros(1,N);
sum_Sxy = zeros(1,N);

for i = (1:K)
	Xd = fft(x(i,:));
	Yd = fft(y(i,:));

	Sxx = 2*dt/N*(Xd.*conj(Xd));
	sum_Sxx = sum_Sxx + fftshift(Sxx);
	Syy = 2*dt/N*(Yd.*conj(Yd));
	sum_Syy = sum_Syy + fftshift(Syy);
	Sxy = 2*dt/N*(Xd.*conj(Yd));
	sum_Sxy = sum_Sxy + fftshift(Sxy);
end

Sxx = sum_Sxx/K;
Syy = sum_Syy/K;
Sxy = sum_Sxy/K;

cohr = Sxy.*conj(Sxy) ./ (Sxx.*Syy);

plot(faxis, cohr)
xlabel('Freq [Hz]'); 
ylabel('Squared Coherence'); 
ylim([-0.1,1.1]);
%
%(Q.iv) Based on your analysis, briefly summarize your results in a few
%sentences.  Specifically address: 1) Do the data exhibit rhythmic activity?
%And 2) Do the signals x and y exhibit coupling?

%% Challenge Problem #3.  Case Study 2.
%Download the data 'Case_Study_2.mat' from Blackobard.  You'll find the
%file contains three variables,
%
%  x = observation of time series x [ntrials, time].
%  y = observation of time series y [ntrials, time].
%  dt = sampling interval in [s].
%
%Both x and y are observed simultaneuously for 100 trials, each trial lasting 1s.
%The sampling interval dt = 0.001s.  Use the tools developed in this lab
%to analyze these data.  Please compute:
%
%(Q.i) The power spectra of x and y.  For these data, compute the power spectra
% for each trial, and average the results over all trials.
%
%(Q.ii) The cross correlation between x and y.  Again, for these data, compute
%the cross correlation for each trial, and average the results over all trials.
%
%(Q.iii) The coherence between x and y.
%
%(Q.iv) Based on your analysis, briefly summarize your results in a few
%sentences.  Specifically address: 1) Do the data exhibit rhythmic activity?
%And 2) Do the signals x and y exhibit coupling?

%% Advanced Challenge Problem #4.
%(Q.i) Examine the MATLAB function 'mscohere' and use it to
%compute the coherence between simulated data.

cohr = zeros(129, 1);
for i = (1:K)
    Cxy = mscohere(x(i,:),y(i,:));
    cohr = cohr + Cxy;
end
cohr = cohr/K;

plot(cohr)
xlabel('Freq [Hz]'); 
ylabel('Squared Coherence'); 
ylim([-0.1,1.1]);
%
%(Q.ii) Examine the function 'coherencyc' in the Chronux Toolbox found here,
%
%    http://www.chronux.org/
%
%and use the function to compute the coherence between *single trial* simulated
%noisy data.
