function y = get_fft(data,varargin)
	fr = dsCalcFR(data,'bin_size',1, 'bin_shift',1);
	fr_soma = {fr.soma_V_FR(1)};
	
	[T_total, numcells] = size(data.soma_V);
    dt = 0.1;
	%numD1s = size(data.D1_V,2);
	%numD2s = size(data.D2_V,2);
    T_total = T_total-1;
    new_T = T_total/10 + 1;
	time = zeros(1,new_T);
    for j = 1:new_T 
        time(j) = (j-1)*dt;
    end
	T_total = size(data.soma_V,1)-1;
    T_start = T_total*0.25; 
	v_new = data.soma_V(T_start:T_total+1,:);

	if numcells > 1
        lfp = mean(v_new');
    else
        lfp = v_new';
    end
    %%%%%%%%%%%%%%%%%%%% spectra
    m = mean(lfp);
    signal = lfp - m; %zero-center
    signal = double(detrend(signal));

	[y,f] =  pmtm(signal,[],[0:150],1000/dt);
    y = {y};
end