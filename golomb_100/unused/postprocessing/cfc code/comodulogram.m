function[h] = comodulogram(d,dt,min_low,max_low,min_hi,max_hi, jnum)
    % Tort et al. Measuring phase-amplitude coupling between neuronal
	% oscillations of different frequencies. JOURNAL OF NEUROPHYSIOLOGY
	% (2010) vol. 104 (2) pp. 1195-210
	
	%input: d (a 1xn vector of data), jnum (the frequency interval)

    N = length(d);  
    Fs = 1/dt;          %Define useful parameters.
    t = (0:N-1)*dt;                     %Define a time axis.
    plot(t, d)
    xlabel('Time [s]')

    periodogram(d, [], [], 1/dt)

    h = zeros(((max_hi - min_hi)/jnum)+1,(max_low-min_low)+1);
    for i = min_low:max_low
        for j = min_hi:jnum:max_hi
            deg=2;                              %Sets the filter order.
            Wn = [i*2/Fs,(i+1)*2/Fs];               %Define the low frequency window of interest.
            [B,A] = butter(deg,Wn,'bandpass');  %Apply the filter to isolate the band.
            dlo = filtfilt(B,A,d);

            Wn = [j*2/Fs,(j+jnum)*2/Fs];            %Define the high frequency window of interest.
            [B,A] = butter(deg,Wn,'bandpass');  %Apply the filter to isolate the band.
            dhi = filtfilt(B,A,d);

            dlo = dlo(N/4:end-N/4-1);
            dhi = dhi(N/4:end-N/4-1);
            taxis = t(N/4:end-N/4-1);

            phi = angle(hilbert(dlo));    %Phase of low frequency signal.
            amp = abs(hilbert(dhi));      %Amplitude envelope of high frequency signal.

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
            h((j/jnum)-(min_hi/jnum) + 1,i-min_low + 1) = max(a_mean)-min(a_mean);                %The difference between the max and min modulation.
			
			%figure(11)
			%plot(p_mean, a_mean, 'k', 'LineWidth', 1);
			%axis tight
			%xlabel('Low frequency phase');  ylabel('High frequency envelope height difference');
			%title(['Metric h=' num2str(h)])
        end
    end

    figure
    imagesc([min_low max_low], [min_hi max_hi], h);
    xlabel('Phase frequency');
    ylabel('Amplitude envelope frequency');
    title(['Comodulogram']);
    colorbar;
end