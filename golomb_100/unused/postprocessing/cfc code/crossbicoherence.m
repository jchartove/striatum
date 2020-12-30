%Compute the bicoherence of data matrix d.
%The format of d is . . . d = d[number of trials, time of trial]
%f0 = The sampling frequency.
%fmax = the maximum frequency of interest.
%
%June 20, 2007.  MAK

function [b, faxis] = crossbicoherence(d1, d2, f0, fmax)

  sz = size(d1);  %Determine the size to define useful parameters.
  K = sz(1);     %Number of trials.
  N = sz(2);     %Number of indices per trial.

  faxis = (0:N/2) / (N)*f0;       %Frequency axis, probably in Hz.
  good = find(faxis < fmax);      %The part of the frequency axis to consider.
  
  b = zeros(length(good), length(good));   %Variable to hold the results.
  faxis = faxis(good);                     %Keep only faxis of interest.
  
  numerator = complex(zeros(length(good)), zeros(length(good)));
  powx = zeros(N,1);
  powy = zeros(N,1);
  
  for k=1:K
      x = fft(hann(N).*d1(k,:)');    %Take the FFT of segment with Hanning.
      y = fft(hann(N).*d2(k,:)');    %Take the FFT of segment with Hanning.
      numTemp = complex(zeros(length(good)), zeros(length(good)));

      for f1=1:length(good)         %For each freq pair, compute numerator.
          for f2=1:length(good)
              numTemp(f1,f2) = x(f1)*x(f2)*conj(y(f1+f2));
          end
      end
      
      numerator = numerator + numTemp / K;   %Compute trial average of numerator.
      powx = powx + x.*conj(x) / K;            %Compute FFT squared, and avg over trials.
      powy = powy + y.*conj(y) / K;            %Compute FFT squared, and avg over trials.
  end

  for f1=1:length(good)                      %Compute bicoherence.
      for f2=1:length(good)
          b(f1,f2) = abs(numerator(f1,f2)) / sqrt(powx(f1)*powx(f2)*powy(f1+f2));
      end
  end

end  