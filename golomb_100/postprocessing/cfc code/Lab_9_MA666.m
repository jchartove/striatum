%%  MA666 - Lab 1 (9):  Correlation and coherence
%   In this lab we will review the power spectrum, develop our own code to
%   compute the coherence and cross correlation, and apply these measures to
%   example simulated data.

%%  Preliminaries.
%   Text preceded by a '%' indicates a 'comment'.  This text should appear
%   green on the screen.  I will use comments to explain what we're doing 
%   and to ask you questions.  Also, comments are useful in your own code
%   to note what you've done (so it makes sense when you return to the code
%   in the future).  It's a good habit to *always* comment your code.  I'll
%   try to set a good example, but won't always . . . 

%%  Part 1:  Warm-up, Computing the power spectrum in 3 ways.
%   We discussed in class the discrete Fourier transform and power
%   spectrum.  Let's compute these quantities in three ways.  We'll begin
%   by writing out every step, and conclude using a built-in MATLAB
%   routine.

%   Let's begin by clearing the workspace.
clear

%  First, we'll create simulated data.
dt = 0.001;             %The sampling interval.
N = 1000;               %The number of time steps.
t = (0:N-1)*dt;         %Define a time axis.
T = N*dt;               %Define the total time of the simulation.

%  Then, create the data, a simple sinusoid + noise.
x = sin(2.0*pi*t*10) + randn(1,N);

%  Now, define the frequencies at which we'll evaluate the Fourier
%  transform.  We consider frequencies between the +/- Nyquist frequency
%  (the sampling frequency divided by two), in steps of 1/T = the frequency
%  resolution.  Here we arrange the frequencies following the MATLAB
%  convention.
fj = [(0:N/2-1)*(1/T), (-N/2:-1)*(1/T)];

%  For each freq compute the Fourier transform and save in "X".
X = zeros(1,length(fj));
n = (1:N);
for j=1:length(fj)
    X(j) = sum(x .* exp(-2*pi*1i*fj(j)*t));
end

%  Compute the power spectrum.
pow = 2*dt^2/T * (X.*conj(X));

%  And plot it on a decibel scale.
subplot(3,1,1)
plot(fj, 10*log10(pow), '*');
xlim([-20 20]);  ylim([-50 0])
xlabel('Freq [Hz]');  ylabel('Power')

%  Second, let's compute the Fourier transform using built-in MATLAB
%  routines.  Then, it's one step,
X = fft(x);

%  Re-compute the power spectrum.
pow = 2*dt^2/T * (X.*conj(X));

%  And plot it on a decibel scale.  Note that we'll use the same frequency
%  axis defined above.
subplot(3,1,2)
plot(fj, 10*log10(pow), '*');
xlim([-20 20]);  ylim([-50 0])
xlabel('Freq [Hz]');  ylabel('Power')

%  Third, use a built-in MATLAB routine to compute and plot spectrum.
subplot(3,1,3)
periodogram(x, [], N, 1/dt);
xlim([-20 20]);  ylim([-50 0])

%  All three measures should give the same result.  Confirm this.

%%  Part 2:  Compute the cross covariance and apply to noisy data.
%   Again, we'll generate artificial data, here starting with
%   two noisy signals.

%   Let's begin by clearing the workspace.
clear

%   And now create the simulated data.
dt = 0.001;                 %The sampling interval.
N = 1000;                   %The number of time steps.
t = (0:N-1)*dt;             %Define a time axis.
x = randn(1,length(t));     %Signal x.
y = randn(1,length(t));     %Signal y.

%   Plot the two signals to have a look.  Visual inspection is always a
%   good place to start data analysis.

subplot(2,1,1)
plot(t,x)
hold on
plot(t,y, 'g')
hold off
xlabel('Time [s]')

%   These signals are *not* coupled (right?  They're just noise).  Let's
%   compute the cross covariance.  To do so, we'll use a built-in MATLAB
%   routine.  Can you already guess what values the cross covariance will
%   obtain?

[rxy, lag] = xcorr(x,y,'biased');

%IN LAB Q:  Note the use of the keyword "biased".  What does the keyword
%do?

%   And plot the results.
subplot(2,1,2)
plot(lag, rxy)
xlabel('Lags [ms]');  ylabel('Cross Covariance');

%%  Part 3:  Relationship of power spectrum to autocovariance.
%   In class, we showed the relationship between the power spectrum of a
%   signal and its autocovariance.  Namely, they are Fourier transform
%   pairs.  In class we showed that,
%
%        The power spectrum is the FT of the auto-covariance.
%
%   It's also true that,
%
%        The auto-covariance is the inverse FT of the power spectrum.
%
%   We'll show this second statement below for simulated data.  

%   Again, let's begin by clearing the workspace.
clear

%   And make a signal we understand,
dt = 0.001;             %The sampling interval.
N = 1000;               %The number of time steps.
t = (0:N-1)*dt;         %Define a time axis.
T = N*dt;               %Define the total time of the simulation.

%   Let's begin with a sinusoid + noise,
x = sin(2*pi*t*10) + randn(1,length(t));

%  Compute the Fourier transform.
X = fft(x);

%  Compute the power spectrum, then take the inverse Fourier transform of
%  the power spectrum, and scale by 1/(2 dt).  The result should be the
%  auto-covariance of x.
pow = 2*dt^2/T* X.*conj(X);
acFFT = 1/(2*dt)*ifft(pow);

%  Rearrange to put negative lags in the correct place.  To do so, we'll
%  use the MATLAB built-in function 'fftshift',
acFFT = fftshift(acFFT);

%  Now, compute the autocovariance using our cross-covariance routine.
[ac, lags] = xcorr(x,x, 'biased');            %Compute covariance.
ac = ac(N/2:end-N/2);                         %Consider lags -N/2 to N/2.
lags = lags(N/2:end-N/2)*dt;                  %Consider lags -N/2 to N/2.

%  To compare, plot the results.
clf()
subplot(2,1,1)
plot(t,x);
ylabel('Data');  xlabel('Time [s]')

subplot(2,1,2)
plot(lags, acFFT);
hold on;
plot(lags, ac, 'go');
hold off
ylabel('Covariance');  xlabel('Lags [s]')

%  Notice we find excellent correspondence between the two methods for
%  computing the autocovariance.  This correspondence is especially good at
%  short time lags (near 0).  You'll notice that the correspondence
%  decreases as the lag progresses.  We can improve this correspondence by
%  performing circular shifts of the time series and computing the
%  auto-covariance, rather than linearly shifting the time series so that
%  points "fall off" the end and are lost to the covaraince calculation.

%%  Part 4:  Compute the the squared coherence (by hand) between two signals.
%   In class we developed an expression for the squared coherence.  We will
%   now write code to compute the squared coherence.  We'll work through
%   the following code in lab, and make sure we understand it.

%   Again, let's begin by clearing the workspace.
clear

%   Let's define two signals.  Can you guess their coherence?
dt = 0.001;             %The sampling interval.
N = 1000;               %The number of time steps.
t = (0:N-1)*dt;         %Define a time axis.
T = N*dt;               %Define the total time of the simulation.
x = randn(1,N);         %Signal x, random noise.
y = randn(1,N);         %Signal y, random noise.

%   Let's also define the frequency axis.  We'll do so using slightly
%   different notation than above.  Note that, here, we've already shifted
%   the frequency axis so that the frequencies increase from negative to
%   positive,
df = 1/T;
fNQ = 1/dt/2;
faxis = (-fNQ:df:fNQ-df);

%  Compute the Fourier transform of x.
Xd = fft(x);

%  Compute the Fourier transform of y.
Yd = fft(y);

%  Compute the auto and cross spectra, and organize vectors by increasing frequency.
Sxx = 2*dt/N*(Xd.*conj(Xd));  Sxx = fftshift(Sxx);
Syy = 2*dt/N*(Yd.*conj(Yd));  Syy = fftshift(Syy);
Sxy = 2*dt/N*(Xd.*conj(Yd));  Sxy = fftshift(Sxy);

%  Compute the squared coherence.
cohr = Sxy.*conj(Sxy) ./ (Sxx.*Syy);

%  Plot the results.
subplot(5,1,1)
plot(t,x)
hold on
plot(t,y, 'g')
hold off
xlabel('Time [s]');  ylabel('Data')

subplot(5,1,2)
plot(faxis, Sxx);  ylabel('Abs(Sxx)^2');  xlim([-50,50])
subplot(5,1,3)
plot(faxis, Syy);  ylabel('Abs(Syy)^2');  xlim([-50,50])
subplot(5,1,4)
plot(faxis, Sxy.*conj(Sxy));  ylabel('Abs(Sxy)^2');  xlim([-50,50])
subplot(5,1,5)
plot(faxis, cohr)
xlabel('Freq [Hz]'); ylabel('Squared Coherence'); xlim([-50,50]); ylim([-0.1,1.1])

% Hmm . . . something odd happening here?
%X_j is a complex number in the complex plane. In polar coordinates: X_j =
%A_j*e^{i*phi_j} where A_j is amplitude, phi_j is phase
%same for Y_j = B_j*e^{i*theta_j}
%Sxx,j = (2*delta^2/T)X_j*X_j^* =
%(2*delta^2/T)(A_j*exp(i*phi_j))(A_j*exp(-i*phi_j))= (2*delta^2/T)A_j^2
%Syy,j = (2*delta^2/T)B_j^2
%Sxy,j = (2*delta^2/T)A_j*B_j*exp(i(phi_j-theta_j))
%K^2_xy,j = (|Sxy,j|^2)/(Sxx,j*Syy,j) =
%(A_j^2*B_j^2*exp(i(phi_j-theta_j))*exp(-i(phi_j-theta_j)))/A_j^2*B_j^2 = 1
%so for all frequencies j the phase difference is 1... this isn't very
%helpful
%to compute coherence we actually need multiple trials - we want to know
%whether x and y are coherent over trials
%K^2_xy,j = (|<Sxy,j>|^2)/(<Sxx,j>*<Syy,j>) where <> is the average across
%trials
%<Sxy,j> = (2*delta^2/T)<X_j^(k)*Y_j^*(k)>_k where k signifies trial number
%= (2*delta^2/T)*something something
%Assume A_j^k = B_j^k = A_j^0. same amplitude for x and y over all trials
%K^2_{xy,j} = |<exp(i*Phi_j)^(k)>_k|^2 = the phase difference over trials
%what is exp(i*Phi_j)^(0)? it is a point on the unit circle on the complex plane
%at frequency j, x and y have a constant phase difference over trials: 
%Phi_j^k = Phi_j^(0) for all k
%if there is strong coherence, there will be a constant phase difference -
%<exp(i*Phi_j)^(0)>_k will be a straight line of length K, divided by K, in
%the complex plane because all of the vectors are in the same direction
%K^2_xy,j = |<exp(i*Phi_j)^(0)>_k|^2 which should be close to 1
%if there is a random phase difference between x and y across 0, then 
%<exp(i*Phi_j)^(0)>_k is all over the unit circle, vectors cancel out, K is
%near 0

