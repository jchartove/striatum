function outdata = gvdelta(data,varargin)
	%y = get_fft(data);
    y_array = dsImportResults(pwd,@get_fft);
    y = cell2mat(y_array(data.simulator_options.sim_id));
	outdata = sum(y(1:3));
end