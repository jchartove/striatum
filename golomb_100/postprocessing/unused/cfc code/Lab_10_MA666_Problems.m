%  Due: Nov 13, 2014
%
%  x.  Please create a PDF with your answers.
%  x.  Please email you answers to:  mak@bu.edu
%  x.  In your email, please use the Subject:  MA666 Lab 2
%  
%  x.  Please answer at least **four** "Challenge Problems".
%  Your solutions should include figures (with axes labeled
%  and captions) and explanatory text.  Please include any
%  MATLAB code you write;  I'll be looking for comments in your
%  code!  The comments should be detailed enough that you could share your
%  code with any of your classmates / colleagues, and he/she could easily
%  read and use your code.
%
%  x.  You may work together.

%% Challenge Problem #1.
%In Part 3 of lab, we considered a single low frequency band and a
%single high frequency band.  We might want to explore a larger number of
%frequency bands.  To do so, we need a *comodulogram*.
%
%Use the code developed in lab to define a new function that computes
%the comodulogram. Your comodulogram should have two axes:
%
%  x-axis:  The phase frequency (e.g., 5Hz to 12Hz in 1Hz steps).
%  y-axis:  The amplitude envelope frequency (e.g., 50Hz to 200Hz in 10Hz steps).
%
%For each pair of (x-axis, y-axis) values, determine "h" (as defined above)
%and plot the (3-dimensional) results.
%
%HINT 1.  See sample comodulograms here,
%         http://www.pnas.org/content/105/51/20517/F5.expansion.html
%HINT 2.  To display h versus your x-axis and y-axis use the MATLAB command,
%         "imagesc"

%% Challenge Problem #2.
%Generate a simulated data set with multiple trials that displays
%strong bicoherence between two frequencies f1 and f2, and their
%sum f1+f2.  HINT:  The bicoherence detects quadratic interactions
%between sinusoids.

%% Challenge Problem #3.
%Download the data set "data_2.mat" from Blackboard.  Analyze this
%signal using all of the tools at your disposal.  Consider:
%
%(Q.i) Visual inspection,
N = length(d);
t = (0:N-1)*dt;
plot(t,d);
fNQ = 1/dt/2;
xlabel('Time (s)'); ylabel('Amplitude');
title('Raw data of data_2.mat');

%(Q.ii) Power spectrum,

T = N*dt; 
fj = [(0:N/2-1)*(1/T), (-N/2:-1)*(1/T)];

powX = zeros(1,N);
X = fft(d);
powX = powX + 2*dt^2/T* X.*conj(X);

plot(fj, powX);
xlim([0 40])
xlabel('Freq [Hz]');  ylabel('Power');
title('Power spectrum of data_2.mat');

plot(fj, 10*log10(powX));
xlabel('Freq [Hz]');  ylabel('Power [dB]');
title('Power spectrum of data_2.mat');

%(Q.iii) Autocorrelation,

[rxx, lag] = xcorr(d,d,'biased');
rhoxx = rxx/(var(d)^2);

plot(lag, rhoxx)
xlabel('Lags [ms]');  ylabel('Autocorrelation');
title('Autocorrelation of data_2.mat');

%(Q.iv) Cross-frequency coupling,
comodulogram('data_2.mat',1,40,40,400)

%,and (Q.v) bicoherence.

[b, faxis] = bicoherence(d, 1/dt, 50);
imagesc(faxis, faxis, b, [0,1])
colorbar
axis xy
xlabel('Frequency [Hz]');  ylabel('Frequency [Hz]');  
title('Bicoherence of data_2.mat');

%For each measure, what do you find?

%%  Challenge Problem #4.
%Download the data 'Case_Study_2.mat' from MA666 Lab 9 (last week's lab),
%and compute:
%
%(Q.i) The bicoherence of x.
%

[b, faxis] = bicoherence(x, 1/dt, 20);
imagesc(faxis, faxis, b, [0,1])
colorbar
axis xy
xlabel('Frequency [Hz]');  ylabel('Frequency [Hz]');  
title('Bicoherence of x');
%(Q.ii) The bicoherence of y.

[b, faxis] = bicoherence(y, 1/dt, 20);
imagesc(faxis, faxis, b, [0,1])
colorbar
axis xy
xlabel('Frequency [Hz]');  ylabel('Frequency [Hz]');  
title('Bicoherence of y');
%
%and (Q.iii) the cross-bicoherence between x and y.  NOTE: You will need to
%update the bicoherence.m function to compute the cross-bicoherence!
%
[b, faxis] = crossbicoherence(x, y, 1/dt, 20);
imagesc(faxis, faxis, b, [0,1])
colorbar
axis xy
xlabel('Frequency [Hz]');  ylabel('Frequency [Hz]');  
title('Cross-bicoherence of x and y');
%(Q.iv) Compare your results to the coherence analysis from Lab 9.  What
%differences do you find?

%%  Advanced Challenge Problem #5.
%Many alternative methods exist to compute CFC.  Choose and implement an
%alternative method from (posted on Blackboard),
%
% Tort et al. Measuring phase-amplitude coupling between neuronal
% oscillations of different frequencies. JOURNAL OF NEUROPHYSIOLOGY
% (2010) vol. 104 (2) pp. 1195-210
%
%When implementing this method, please apply it to simulated data.  You may
%simulate your own data, or use "data_1.mat" from Lab 10.

%%  Advanced Challenge Problem #6.
%Another alternative method to compute the CFC is presented in (posted on
%Blackboard),
%
% Kramer and Eden. Assessment of cross-frequency coupling with confidence
% using generalized linear models. J Neurosci Methods (2013) vol. 220 (1)
% pp. 64-74
%
%Using the MATLAB code provided in this article (posted on Blackboard as
%the file "GLM_CFC_for_paper.m" and see Appendix A of the article), compute
%the CFC of "data_1.mat" from Lab 10 using the proposed method.  Explain
%your results!

%% Challenge Problem #7.
%Make up your own challenge.  Make it interesting.

